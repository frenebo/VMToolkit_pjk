#ifndef __NEW_FORCE_EFIELD_ON_CELL_BOUNDARY_PIXELATED_HPP__
#define __NEW_FORCE_EFIELD_ON_CELL_BOUNDARY_PIXELATED_HPP__

#include <optional>
#include "force.hpp"



namespace VMSim
{
namespace PixelatedElectricStuff
{
    using std::optional;
    using std::cout;
    using std::endl;
    
    class PixelIntersectionsDescriptor
    {
    public:
        PixelIntersectionsDescriptor()
        {
            // aaa
        }
        
        bool has_top_intersection() const { return _top_intersection_loc.has_value(); }
        bool has_bottom_intersection() const { return _bottom_intersection_loc.has_value(); }
        bool has_left_intersection() const { return _left_intersection_loc.has_value(); }
        bool has_right_intersection() const { return _right_intersection_loc.has_value(); }
        
        double top_intersection_loc() const { return _top_intersection_loc.value(); }
        double bottom_intersection_loc() const { return _bottom_intersection_loc.value(); }
        double left_intersection_loc() const { return _left_intersection_loc.value(); }
        double right_intersection_loc() const { return _right_intersection_loc.value(); }
        
        void set_top_intersection_loc(double loc) {
            if (_top_intersection_loc.has_value()) { throw runtime_error("_top_side_intersection_loc already has value"); }
            _top_intersection_loc = loc;
        }
        
        void set_bottom_intersection_loc(double loc) {
            if (_bottom_intersection_loc.has_value()) { throw runtime_error("_bottom_side_intersection_loc already has value"); }
            _bottom_intersection_loc = loc;
        }
        
        void set_left_intersection_loc(double loc) {
            if (_left_intersection_loc.has_value()) { throw runtime_error("_left_side_intersection_loc already has value"); }
            _left_intersection_loc = loc;
        }
        
        void set_right_intersection_loc(double loc) {
            if (_right_intersection_loc.has_value()) { throw runtime_error("_right_side_intersection_loc already has value"); }
            _right_intersection_loc = loc;
        }
        
        int count_intersections() const {
            int n_intersects = 0;
            
            if (has_top_intersection()) n_intersects++;
            if (has_bottom_intersection()) n_intersects++;
            if (has_left_intersection()) n_intersects++;
            if (has_right_intersection()) n_intersects++;
            
            return n_intersects;
        }
    private:
        optional<double> _top_intersection_loc;
        optional<double> _bottom_intersection_loc;
        optional<double> _left_intersection_loc;
        optional<double> _right_intersection_loc;
    };
    
    class GridCoord
    {
    public:
        GridCoord(int x, int y): _x{x}, _y{y}
        {
            // pass
        }
        GridCoord(): _x{0}, _y{0}
        {
            // pass
        }
        
		//! Test equality
		bool operator==(const GridCoord &gc)
		{
			return (_x == gc._x && _y == gc._y);
		}
        
        bool operator!=(const GridCoord& gc)
        {
            return ! operator==(gc);
        }
        
        int x() const { return _x; }
        int y() const { return _y; }
    private:
        int _x;
        int _y;
    };
    
    class GridSpec
    {
    public:
        GridSpec(
            double origin_x,
            double origin_y,
            double spacing_x,
            double spacing_y,
            int ncells_x,
            int ncells_y
        ):  _origin_x{origin_x},
            _origin_y{origin_y},
            _spacing_x{spacing_x},
            _spacing_y{spacing_y},
            _ncells_x{ncells_x},
            _ncells_y{ncells_y}
        {
            // nothing
        }
        
        double origin_x() const { return _origin_x; }
        double origin_y() const { return _origin_y; }
        double spacing_x() const { return _spacing_x; }
        double spacing_y() const { return _spacing_y; }
        int ncells_x() const { return _ncells_x; }
        int ncells_y() const { return _ncells_y; }
    private:
        double _origin_x;
        double _origin_y;
        double _spacing_x;
        double _spacing_y;
        int _ncells_x;
        int _ncells_y;
    };
    
      
    /*
        * Represents the parameters for a given face (cell) - e.g. its electric charge.
        * This will define how it interacts and experiences forces from the field.
    */
    class EFieldPixFaceConfig
    {
    public:
        EFieldPixFaceConfig(
            double charge
        ): _charge{charge}
        {
            // empty 
        }
        
        // default constructor
        EFieldPixFaceConfig(): _charge{0}
        {
        }
        
        double charge() const { return _charge; }
    private:
        double _charge;
    };
    
    class EdgeComputationResult
    {
    public:
        EdgeComputationResult(
            Vec integrated_efield_with_respect_to_dlength
        ): _int_efield_wrt_dl{integrated_efield_with_respect_to_dlength}
        {
            // nothing
        }
        
        EdgeComputationResult(): _int_efield_wrt_dl{0.0,0.0}
        {
            // nothing
        }
        
        Vec integrated_efield_with_respect_to_dlength() const { return _int_efield_wrt_dl; }
    private:
        Vec _int_efield_wrt_dl;
    };
    
    
    class NEWForceEFieldOnCellBoundPixelated : public Force
    {
    public:
        // The actual class definition
		NEWForceEFieldOnCellBoundPixelated(const System &sys) : Force{sys}
        {
		}
        
		virtual ~NEWForceEFieldOnCellBoundPixelated() {}
        
        
        void set_global_params(
            const params_type& num_params,
            const map<string,string>& str_params,
            const std::map<string, int>& int_params,
            const map<string, vector<double>> flt_array_params,
            bool verbose
        ) override
        {
            if (verbose) {
                cout << "  NEWForceEFieldOnCellBoundPixelated::set_global_params - starting" << endl;
            }
            _check_for_param_names(
                num_params,
                str_params,
                int_params,
                flt_array_params,
                {
                    "grid_origin_x",
                    "grid_origin_y",
                    "grid_spacing_x",
                    "grid_spacing_y"
                },
                {},
                {
                    "grid_ncells_x",
                    "grid_ncells_y"
                },
                {
                    "flattened_field_vals_by_pixel"
                }
            );
            
            _gridspec = GridSpec(
                num_params.at("grid_origin_x"),
                num_params.at("grid_origin_y"),
                num_params.at("grid_spacing_x"),
                num_params.at("grid_spacing_y"),
                int_params.at("grid_ncells_x"),
                int_params.at("grid_ncells_y")
            );
            
            // The values of the electric field grid
            const size_t ncells_x = _gridspec.value().ncells_x();
            const size_t ncells_y = _gridspec.value().ncells_y();
            
            const vector<double>& flattened_field_vals_elements = flt_array_params.at("flattened_field_vals_by_pixel");
            size_t tot_num_pixels = ncells_x * ncells_y;
            if (flattened_field_vals_elements.size() != tot_num_pixels * 2) {
                cout << "flattened_field_vals_elements.size() - " << flattened_field_vals_elements.size() << endl;
                cout << "_gridspec.ncells_x() - " << ncells_x << endl;
                cout << "_gridspec.ncells_y() - " << ncells_y << endl;
                
                throw runtime_error("Invalid size for 'flattened_field_vals_by_pixel' - size does not match ncells_x * ncells_y * 2 (2 for x and y)");
            }
            
            _flattened_field_vecs = vector<Vec>();
            _flattened_field_vecs.value().resize(tot_num_pixels, Vec(0.0,0.0));
            
            if (verbose) {
                cout << "Populating the flattened field vector" << endl;
            }
            // Construct rows
            for (size_t flat_pos=0; flat_pos < tot_num_pixels; flat_pos++) {
                _flattened_field_vecs.value().at(flat_pos) = Vec(
                    flattened_field_vals_elements.at(flat_pos*2),
                    flattened_field_vals_elements.at(flat_pos*2 + 1)
                );
            }
        }
        
        Vec get_field_at_gridpos(size_t col, size_t row, bool verbose) const
        {
            if (!_flattened_field_vecs) {
                throw runtime_error("get_field_at_pos - Cannot get field - field data has not been set yet");
            }
            if (!_gridspec) {
                throw runtime_error("get_field_at_pos - cannot get field, gridspec has not been set yet");
            }
            
            if (col >= _gridspec->ncells_x()) {
                cout << "col, row:" << col << "," << row << endl;
                throw runtime_error("Invalid column - exceeds bounds");
            }
            if (row >= _gridspec->ncells_y()) {
                cout << "col, row:" << col << "," << row << endl;
                throw runtime_error("Invalid row - exceeds bounds");
            }
            size_t flat_pos = col * _gridspec->ncells_y() + row;
            if (flat_pos > _flattened_field_vecs->size()) {
                throw runtime_error("??? this shouldn't happen");
            }
            
            return _flattened_field_vecs->at(flat_pos);
        }
        
        void set_face_params_facewise(
            const vector<int>& fids,
            const vector<params_type>& params,
            bool verbose
        ) override
        {
            _face_params.clear();
            if (fids.size() != params.size()) {
                throw runtime_error("Size of fids and params don't match");
            }
            
            for (int i = 0; i < fids.size(); i++) {
                int fid = fids.at(i);
                const params_type& p = params.at(i);
                
                double cell_charge = p.at("charge");
                EFieldPixFaceConfig cell_params = EFieldPixFaceConfig(cell_charge);
                
                _face_params[fid] = cell_params;
            }
        }
        
		void compute_all_vertex_forces(vector<Vec>& res, bool verbose) override;
        
    private:
        void _check_for_param_names(
            const params_type& num_params,
            const map<string,string>& str_params,
            const std::map<string, int>& int_params,
            const map<string, vector<double>> flt_array_params,
            
            const vector<string>& expected_num_params,
            const vector<string>& expected_str_params,
            const vector<string>& expected_int_params,
            const vector<string>& expected_flt_array_params
        ) const {
            for (const string& name : expected_num_params) {
                if (! num_params.contains(name))
                    throw runtime_error("Missing num param: " + name);
            }
            for (const string& name : expected_str_params) {
                if (! str_params.contains(name))
                    throw runtime_error("Missing str param: " + name);
            }
            for (const string& name : expected_int_params) {
                if (! int_params.contains(name))
                    throw runtime_error("Missing int param: " + name);
            }
            for (const string& name : expected_flt_array_params) {
                if (! flt_array_params.contains(name))
                    throw runtime_error("Missing flt_array param: " + name);
            }
        }
        
        void _clear_compute_cache(bool verbose);
        
        void _cache_computations(bool verbose);
        
        EdgeComputationResult _integrate_field_over_edge(const Edge& edge, bool verbose) const;
        
        vector<GridCoord> _get_edge_pixel_intersections(const Edge& edge, bool verbose) const;
        vector<double> _get_edge_lengths_passing_thru_pixels(const Edge& edge, const vector<GridCoord>& pixels_intersected_by_edge, bool verbose) const;
        Vec _get_vec_coords_for_grid_pos(GridCoord gc, bool verbose) const;
        
        GridCoord _get_grid_coords_for_vector(Vec vec, bool verbose) const;
        
        vector<int> _get_crossings_generalized(int start_coord, int end_coord, bool verbose) const;
        vector<int> _get_column_crossings_of_edge(GridCoord edge_start, GridCoord edge_end, bool verbose) const;
        vector<int> _get_row_crossings_of_edge(GridCoord edge_start, GridCoord edge_end, bool verbose) const;
        double get_relative_position_of_crossing_along_edge(int crossing_coord, double edge_start_coord, double edge_end_coord) const;
        
        double _get_intersection_with_row_and_find_rel_pos_in_column(const Vec& edge_start_VEC, const Vec& edge_end_VEC, int row_to_intersect, int snap_column, bool verbose) const;
        
        double _get_intersection_with_column_and_find_rel_pos_in_row(const Vec& edge_start_VEC, const Vec& edge_end_VEC, int column_to_intersect, int snap_row, bool verbose) const;
        
        vector<Vec> _get_pixel_intersection_absolute_coordinates(const GridCoord& pixel_gridpos, const PixelIntersectionsDescriptor& pix_intersections, bool verbose) const;
        
        optional<GridSpec> _gridspec;
        optional<vector<Vec>> _flattened_field_vecs;
        map<int, EFieldPixFaceConfig> _face_params;
        map<int, EdgeComputationResult> _cached_edge_computation_results;
    };
}
}

#endif
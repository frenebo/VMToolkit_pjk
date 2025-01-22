#ifndef __OLD_FORCE_EFIELD_ON_CELL_BOUNDARY_PIXELATED_HPP__
#define __OLD_FORCE_EFIELD_ON_CELL_BOUNDARY_PIXELATED_HPP__

#include <optional>

#include "force.hpp"

using std::map;
using std::optional;
using std::cout;
using std::map;
using std::endl;
using std::runtime_error;
using std::vector;
using std::string;


namespace VMSim
{
    
    
    
    
	// Force on a vertex
	class OLDForceEFieldOnCellBoundPixelated : public Force
	{
        // The child classes used to represent various data:
	public:
        class PixIntersectionsSpec {
        public:
            PixIntersectionsSpec()
            {
            }
            
            void add_left_intersection(float ypos) {
                if (_left_intersection_pos) { throw runtime_error("Already set _left_intersection_pos"); }
                _left_intersection_pos = ypos;
            }
            
            void add_right_intersection(float ypos) {
                if (_right_intersection_pos) { throw runtime_error("Already set _right_intersection_pos"); }
                _right_intersection_pos = ypos;
            }
            
            void add_top_intersection(float pos) {
                if (_top_intersection_pos) { throw runtime_error("Already set _top_intersection_pos"); }
                _top_intersection_pos = pos;
            }
            
            void add_bot_intersection(float pos) {
                if (_bot_intersection_pos) { throw runtime_error("Already set _bot_intersection_pos"); }
                _bot_intersection_pos = pos;
            }
            
            int count_intersections() const {
                int cnt = 0;
                
                if (_left_intersection_pos) cnt++;
                if (_right_intersection_pos) cnt++;
                if (_bot_intersection_pos) cnt++;
                if (_top_intersection_pos) cnt++;
                
                return cnt;
            }
            
            const optional<float>& left_intersection() const { return _left_intersection_pos; }
            const optional<float>& right_intersection() const { return _right_intersection_pos; }
            const optional<float>& top_intersection() const { return _top_intersection_pos; }
            const optional<float>& bot_intersection() const { return _bot_intersection_pos; }
        private:
            optional<float> _left_intersection_pos;
            optional<float> _right_intersection_pos;
            optional<float> _top_intersection_pos;
            optional<float> _bot_intersection_pos;
        };
        /*
            * A container to represent the grid of the force field
        */
        class GridSpec {
        public:
            GridSpec(
                float orig_x,
                float orig_y,
                float g_space,
                int nx,
                int ny
            ):
                _origin_x{orig_x},
                _origin_y{orig_y},
                _grid_spacing{g_space},
                _ncells_x{nx},
                _ncells_y{ny}
            {
                //
            }
            float origin_x() const {return _origin_x; }
            float origin_y() const {return _origin_y; }
            float grid_spacing() const {return _grid_spacing; }
            int ncells_x() const {return _ncells_x; }
            int ncells_y() const {return _ncells_y; }
        private:
            float _origin_x;
            float _origin_y;
            float _grid_spacing;
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
    
    public:
        // The actual class definition
		OLDForceEFieldOnCellBoundPixelated(const System &sys) : Force{sys}
        {
		}
        
		virtual ~OLDForceEFieldOnCellBoundPixelated() {}
      
		void set_global_params(
            const params_type& num_params,
            const map<string,string>& str_params,
            const std::map<string, int>& int_params,
            const map<string, vector<double>> flt_array_params,
            bool verbose
        ) override
        {
            _gridspec = GridSpec(
                num_params.at("grid_origin_x"),
                num_params.at("grid_origin_y"),
                num_params.at("grid_spacing"),
                int_params.at("grid_ncells_x"),
                int_params.at("grid_ncells_y")
            );
            
            // The values of the electric field grid
            const size_t ncells_x = _gridspec.value().ncells_x();
            const size_t ncells_y = _gridspec.value().ncells_y();
            
            const vector<double>& flattened_field_vals_pixels = flt_array_params.at("flattened_field_vals_by_pixel");
            if (flattened_field_vals_pixels.size() != ncells_x * ncells_y * 2) {
                cout << "flattened_field_vals_pixels.size() - " << flattened_field_vals_pixels.size() << endl;
                cout << "_gridspec.ncells_x() - " << ncells_x << endl;
                cout << "_gridspec.ncells_y() - " << ncells_y << endl;
                
                throw runtime_error("Invalid size for 'flattened_field_vals_by_pixel' - size does not match ncells_x * ncells_y * 2 (2 for x and y)");
            }
            
            _field_vals_2d = vector<vector<Vec>>();
            
            // Construct rows
            for (size_t col_idx = 0; col_idx < ncells_x; col_idx++) {
                vector<Vec> field_vals_column;
                field_vals_column.resize(ncells_y, Vec(0.0,0.0));
                
                size_t col_position_in_flattened = col_idx * ncells_y * 2;
                
                for (size_t row_idx = 0; row_idx < ncells_y; row_idx++) {
                    size_t pix_position_in_flattened =  col_position_in_flattened + row_idx * 2;
                    
                    double field_xval = flattened_field_vals_pixels.at(pix_position_in_flattened);
                    double field_yval = flattened_field_vals_pixels.at(pix_position_in_flattened + 1);
                    field_vals_column.at(row_idx) = Vec(
                        field_xval,
                        field_yval
                    );
                    
                    cout << "X: " << field_xval << ", Y: " << field_yval << endl;
                }
                
                _field_vals_2d.value().push_back(field_vals_column);
            }
        }
        
		
		void compute_all_vertex_forces(vector<Vec>& res, bool verbose) override;

		void set_face_params_facewise(const vector<int>& fids, const vector<params_type>& params, bool verbose) override
        {
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
        
    private:
        // void _dump_and_error
        void _clear_compute_cache(bool verbose);
        
        void _cache_all_computations(bool verbose);
        
        Vec _get_pixel_origin_loc(std::pair<int,int> pix_gridpos) const;
        
        std::pair<int, int> _get_pixel_indices(const Vec& loc, bool verbose) const;
        
        double _lineseg_solve_for_y(double x, const Vec& line_start, const Vec& line_end, bool verbose) const;
        
        double _lineseg_solve_for_x(double y, const Vec& line_start, const Vec& line_end, bool verbose) const;
        
        vector<Vec> _get_pixspec_relative_intersection_positions(const PixIntersectionsSpec& pixspec, bool verbose) const;
        
        double _get_intersection_of_line_passthrough_pixel(const PixIntersectionsSpec& pixspec, bool verbose) const;
        
        double _get_intersection_of_line_ending_inside_pixel(
            std::pair<int,int> pix_gridpos,
            const PixIntersectionsSpec& pixspec, 
            const Edge& edge,
            bool verbose
        ) const;
        
        Vec _integrate_field_over_edge_wrt_length(const Edge& edge, bool verbose) const;

        std::map<std::pair<int,int>, PixIntersectionsSpec> _get_edge_pixel_intersections(
            const Edge& edge,
            bool verbose
        ) const;
        
        
        std::optional<GridSpec> _gridspec;
        std::optional< vector<vector<Vec>> > _field_vals_2d;
        std::map<int, EFieldPixFaceConfig> _face_params;
        
        std::map<int, Vec> _cached_edge_integrated_field_wrt_length;
	};
    
}



#endif
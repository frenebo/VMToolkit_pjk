#ifndef __FORCE_EFIELD_ON_CELL_BOUNDARY_PIXELATED_HPP__
#define __FORCE_EFIELD_ON_CELL_BOUNDARY_PIXELATED_HPP__

#include <optional>

#include "force.hpp"

using std::map;
using std::optional;


namespace VMSim
{
	// Force on a vertex
	class ForceEFieldOnCellBoundPixelated : public Force
	{
        // class CellInt
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
                float nx,
                float ny
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
            float ncells_x() const {return _ncells_x; }
            float ncells_y() const {return _ncells_y; }
        private:
            float _origin_x;
            float _origin_y;
            float _grid_spacing;
            float _ncells_x;
            float _ncells_y;
        };
        
        /*
         * Represents the parameters for a given face (cell) - e.g. its electric charge.
         * This will define how it interacts and experiences forces from the field.
        */
        class EFieldFaceParams
        {
        public:
            EFieldFaceParams(
                double charge
            ): _charge{charge}
            {
                // empty 
            }
            
            // default constructor
            EFieldFaceParams(): _charge{0}
            {
            }
            
            double charge() const { return _charge; }
        private:
            double _charge;
        };

        
	public:
		ForceEFieldOnCellBoundPixelated(const System &sys) : Force{sys}
        {
		}
        
		virtual ~ForceEFieldOnCellBoundPixelated() {}
      
		void set_global_params(const params_type& num_params, const map<string,string>& str_params, bool verbose) override
        {
            _gridspec = GridSpec(
                num_params.at("grid_origin_x"),
                num_params.at("grid_origin_y"),
                num_params.at("grid_spacing"),
                num_params.at("grid_ncells_x"),
                num_params.at("grid_ncells_y")
            );
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
                EFieldFaceParams cell_params = EFieldFaceParams(cell_charge);
                
                _face_params[fid] = cell_params;
            }
        }
    private:
		Vec _compute_he_force_for_face(
			const Vertex& vert,
			const HalfEdge& he,
			const Face& f,
			// double edge_len_within_polygon,
			bool verbose=false
		);
        
        Vec _compute_he_force_for_face(
            const Vertex& vert,
            const HalfEdge& he,
            const Face& f,
            bool verbose=false
        ) const;
        
        void _clear_compute_cache(bool verbose);
        void _cache_all_computations(bool verbose);
        void _cache_comp_electric_field_on_edges(bool verbose);
        Vec _get_pixel_origin_loc(std::pair<int,int> pix_gridpos) const;
        std::pair<int, int> _get_pixel_indices(Vec loc) const;
        void _lineseg_solve_for_y(double x, const Vec& line_start, const Vec& line_end) const;
        void _lineseg_solve_for_x(double y, const Vec& line_start, const Vec& line_end) const;
        vector<Vec> _get_pixspec_relative_intersection_positions(const PixIntersectionsSpec& pixspec) const;
        
        double _get_intersection_of_line_passthrough_pixel(const PixIntersectionsSpec& pixspec) const;
        double _get_intersection_of_line_ending_inside_pixel(
            std::pair<int,int> pix_gridpos,
            const PixIntersectionsSpec& pixspec, 
            const Edge& edge
        ) const;
        Vec _integrate_field_over_edge(const Edge& edge, bool verbose) const;

        std::map<std::pair<int,int>, PixIntersectionsSpec> _get_edge_pixel_intersections(
            const Edge& edge,
            bool verbose
        ) const;
        
        
        std::optional<GridSpec> _gridspec;
        std::optional< vector<vector<Vec>> > _field_vals_2d;
        std::map<int, EFieldFaceParams> _face_params;
	};
    
}



#endif
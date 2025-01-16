#ifndef __FORCE_EFIELD_ON_CELLLBOUND_PIXELATED_HPP__
#define __FORCE_EFIELD_ON_CELLLBOUND_PIXELATED_HPP__

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
                float charge
            ): _charge{charge}
            {
                // empty 
            }
            
            float charge() const { return _charge; }
        private:
            float _charge;
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
                
                _face_params[fid] = EFieldFaceParams(
                    p.at("charge"),
                );
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
        
		void _clear_compute_cache(bool verbose);
        void _cache_all_computations(bool verbose);
        
        optional<GridSpec> _gridspec;
        optional< vector<vector<Vec>> > _field_vals_2d;
        map<int, EFieldFaceParams> _face_params;
		// bool _vec_outside_polygon_bbox(const Vec& v, bool verbose);
        // bool _enabled_for_faceidx(int fid, bool verbose);
		// Vec _compute_he_force(const Vertex& vert, const HalfEdge& he, bool verbose=false);
		
        
        
        
		// void _cache_mesh_computations(bool verbose);
		// double _lazy_load_cell_edge_intersection(int vid_from, int vid_to, bool verbose);
		
		
		// PolygonZone _poly_zone;
		// double _E_x_param;
		// double _E_y_param;
		
        // vector<bool> _force_enabled_mask_by_face_index;
		// vector<double> _cell_charges_by_face_index;
		
		// vector<bool> _cached_vertices_outside_bbox;
		// vector<bool> _cached_enabled_for_fid;
		// vector<double> _cached_face_perims;
		// map<std::pair<int, int>, double> _cached_edgelength_intersecting_polygon_by_vtx_ids;
	};
    
}



#endif
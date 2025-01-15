#ifndef __FORCE_EFIELD_ON_CELLLBOUND_PIXELATED_HPP__
#define __FORCE_EFIELD_ON_CELLLBOUND_PIXELATED_HPP__

#include <optional>

#include "force.hpp"


namespace VMSim
{
    class GridSpec {
    public:
        GridSpec(
            float orig_x,
            float orig_y,
            float g_space,
            float nx,
            float ny
        ):
            origin_x{orig_x},
            origin_y{orig_y},
            grid_spacing{g_space},
            ncells_x{nx},
            ncells_y{ny}
        {
            //
        }
        
        float origin_x;
        float origin_y;
        float grid_spacing;
        float ncells_x;
        float ncells_y;
    };

	// Force on a vertex
	class ForceEFieldOnCellBoundPixelated : public Force
	{
        
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
            
        }
    private:
		Vec _compute_he_force_for_face(
			const Vertex& vert,
			const HalfEdge& he,
			const Face& f,
			double edge_len_within_polygon,
			bool verbose=false
		);
        
        std::optional<GridSpec> _gridspec;
		// bool _vec_outside_polygon_bbox(const Vec& v, bool verbose);
        // bool _enabled_for_faceidx(int fid, bool verbose);
		// Vec _compute_he_force(const Vertex& vert, const HalfEdge& he, bool verbose=false);
		
		// void _clear_compute_cache(bool verbose);
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
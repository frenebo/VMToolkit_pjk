/*!
 * \file force_const_vertex_propulsion.hpp
 * \author Paul Kreymborg, frenebo@gmail.com
 * \date 20-Dec-2024
 * \brief ForceEFieldOnCellBoundary class 
*/ 

#ifndef __FORCE_E_FIELD_ON_CELL_BOUNDARY_UNIFORM_HPP__
#define __FORCE_E_FIELD_ON_CELL_BOUNDARY_UNIFORM_HPP__

#include "force.hpp"
#include "polygon_zone.hpp"

#include <utility> // std::pair
#include <string>

namespace VMSim
{

	// Force on a vertex
	class ForceEFieldOnCellBoundary : public Force
	{
	public:
		ForceEFieldOnCellBoundary(const System &sys) : Force{sys}, _E_x_param{0},_E_y_param{0} {
		}
        
		virtual ~ForceEFieldOnCellBoundary() {}
      
		void set_global_params(
			const params_type& num_params,
			const std::map<string,string>& str_params,
			const std::map<string, int>& int_params,
			const map<string, vector<double>> flt_array_params,
			bool verbose
		) override;
		
		void compute_all_vertex_forces(vector<Vec>& res, bool verbose) override;

		void set_face_params_facewise(const vector<int>& fids, const vector<params_type>& params, bool verbose) override;    
    private:
		Vec _compute_he_force_for_face(
			const Vertex& vert,
			const HalfEdge& he,
			const Face& f,
			double edge_len_within_polygon,
			bool verbose=false
		);
		
		bool _vec_outside_polygon_bbox(const Vec& v, bool verbose);
        bool _enabled_for_faceidx(int fid, bool verbose);
		Vec _compute_he_force(const Vertex& vert, const HalfEdge& he, bool verbose=false);
		
		void _clear_compute_cache(bool verbose);
		void _cache_mesh_computations(bool verbose);
		double _lazy_load_cell_edge_intersection(int vid_from, int vid_to, bool verbose);
		
		
		PolygonZone _poly_zone;
		double _E_x_param;
		double _E_y_param;
		
        vector<bool> _force_enabled_mask_by_face_index;
		vector<double> _cell_charges_by_face_index;
		
		vector<bool> _cached_vertices_outside_bbox;
		vector<bool> _cached_enabled_for_fid;
		vector<double> _cached_face_perims;
		map<std::pair<int, int>, double> _cached_edgelength_intersecting_polygon_by_vtx_ids;
	};

}

#endif

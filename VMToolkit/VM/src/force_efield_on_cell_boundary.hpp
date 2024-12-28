/*!
 * \file force_const_vertex_propulsion.hpp
 * \author Paul Kreymborg, frenebo@gmail.com
 * \date 20-Dec-2024
 * \brief ForceEFieldOnCellBoundary class 
*/ 

#ifndef __FORCE_E_FIELD_ON_CELL_BOUNDARY_HPP__
#define __FORCE_E_FIELD_ON_CELL_BOUNDARY_HPP__

#include "force.hpp"
#include "polygon_zone.hpp"
#include <string>

namespace VMTutorial
{

	// Force on a vertex
	class ForceEFieldOnCellBoundary : public Force
	{
	public:
		ForceEFieldOnCellBoundary(const System &sys) : Force{sys}, _E_x{0},_E_y{0} {
		}
        
		virtual ~ForceEFieldOnCellBoundary() {}

		// // computes force on vertex by given half edge
		// Vec compute_he_force(const Vertex<Property> &, const HalfEdge<Property> &, bool verbose=false) override;

      
		void set_global_params(const params_type& num_params, const map<string,string>& str_params, bool verbose) override;
		
		
		Vec compute_he_force_for_face(
			const Vertex<Property>& vert,
			const HalfEdge<Property>& he,
			const Face<Property>& f,
			double edge_len_within_polygon,
			bool verbose=false
		);
		
		bool vec_outside_polygon_bbox(const Vec& v);
		
		Vec compute_he_force(const Vertex<Property>& vert, const HalfEdge<Property>& he, bool verbose=false) override;

		void set_face_params_facewise(const vector<int>& fids, const vector<params_type>& params, bool verbose) override;
    
        bool enabled_for_faceidx(int fid, bool verbose);
    
    private:
        // vector<Vec> _force_by_vidx;
		// vector<Vec> _polygon_vertices;
		PolygonZone _poly_zone;
		double _E_x;
		double _E_y;
		
        vector<bool> _force_enabled_mask_by_face_index;
		vector<double> _cell_charges_by_face_index;
	};

}

#endif

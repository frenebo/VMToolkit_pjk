/*!
 * \file force_const_vertex_propulsion.hpp
 * \author Paul Kreymborg, frenebo@gmail.com
 * \date 17-Dec-2024
 * \brief ForceConstVertexPropulsion class 
*/ 

#ifndef __FORCE_CONST_VERTEX_PROP_HPP__
#define __FORCE_CONST_VERTEX_PROP_HPP__

#include "force.hpp"

namespace VMTutorial
{

	// Force on a vertex
	class ForceConstVertexPropulsion : public Force
	{
	public:
		ForceConstVertexPropulsion(const System &sys) : Force{sys} {
		}
        
		virtual ~ForceConstVertexPropulsion() {}

		// computes force on vertex by a given edge
		Vec compute_he_force(const Vertex<Property> &, const HalfEdge<Property> &, bool verbose=false) override;

		void set_vertex_params_vertexwise(const vector<int>& vids, const vector<params_type>& params, bool verbose) override;
    
        bool enabled_for_vertexidx(int vid, bool verbose);
    
    private:
        vector<Vec> _force_by_vidx;
        vector<bool> _force_enabled_mask_by_vertex_index;
	};

}

#endif

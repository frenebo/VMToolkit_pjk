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
		
		void compute_all_vertex_forces(vector<Vec>& res, bool verbose) override;

		void set_vertex_params_vertexwise(const vector<int>& vids, const vector<params_type>& params, bool verbose) override;
    private:
		Vec compute_he_force(const Vertex &, const HalfEdge &, bool verbose=false);
        bool enabled_for_vertexidx(int vid, bool verbose);
        
		vector<Vec> _force_by_vidx;
        vector<bool> _force_enabled_mask_by_vertex_index;
	};

}

#endif

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
		ForceConstVertexPropulsion(System &sys) : Force{sys} {
		}
        
		virtual ~ForceConstVertexPropulsion() {}

		// computes force on vertex by a given edge
		Vec compute_he_force(const Vertex<Property> &, const HalfEdge<Property> &, bool verbose=false) override;

		// double tension(const HalfEdge<Property> &) override;

		// // Energy calculation
		// double energy(const Face<Property> &) override;
		

		void set_vertex_params_vertexwise(const vector<int>& vids, const vector<params_type>& params, bool verbose) override {
			// Expand the params if necessary
			int max_vid = 0;
			for (const auto& vid : vids) {
                if (vid > max_vid) max_vid = vid;
			}
			
			if (max_vid >= _force_enabled_mask_by_vertex_index.size()) {
				_force_enabled_mask_by_vertex_index.resize(max_vid + 1, false);
				
                // Vec blank_vec = ;
                _force_by_vidx.resize(max_vid + 1, Vec(0.0,0.0));
				
			}
			
			// Set the param values in the vectors
			for (size_t i = 0; i < vids.size(); ++i) {
				int vid = vids[i];
				const params_type& vparam = params.at(i);
                
                _force_enabled_mask_by_vertex_index.at(vid) = true;
                _force_by_vidx.at(vid) = Vec(
                    vparam.at("f_x"),
                    vparam.at("f_y")
                );
			}
		}
    
        bool enabled_for_vertexidx(int vid, bool verbose) {
            if (vid >= _force_enabled_mask_by_vertex_index.size()) return false;
            
            bool is_enabled = _force_enabled_mask_by_vertex_index[vid];
            
            if (verbose && is_enabled) {
                cout << "            Confirmed that constant propulsion force is on for vid=" << vid << endl;
            }
            
            return is_enabled;
        }
    
    private:
        vector<Vec> _force_by_vidx;
        vector<bool> _force_enabled_mask_by_vertex_index;
	};

}

#endif

/*!
 * \file force_const_vertex_propulsion.hpp
 * \author Paul Kreymborg, frenebo@gmail.com
 * \date 20-Dec-2024
 * \brief ForceEFieldOnCellBoundary class 
*/ 

#ifndef __FORCE_E_FIELD_ON_CELL_BOUNDARY_HPP__
#define __FORCE_E_FIELD_ON_CELL_BOUNDARY_HPP__

#include "force.hpp"

namespace VMTutorial
{

	// Force on a vertex
	class ForceEFieldOnCellBoundary : public Force
	{
	public:
		ForceEFieldOnCellBoundary(System &sys) : Force{sys}, _E_x{0},_E_y{0} {
		}
        
		virtual ~ForceEFieldOnCellBoundary() {}

		// computes force on vertex by given half edge
		Vec compute_he_force(const Vertex<Property> &, const HalfEdge<Property> &, bool verbose=false) override;

      
		void set_global_params(const params_type& params, bool verbose) override
		{
			// params.at("n_")/
			string field_type = params.at("field_type");
			
			float E_x = params.at("E_x");
			float E_y = params.at("E_y");
			
			// int num_vertices
			string region_type = params.at("region_type")
			
			// int 
			
			_E_y = E_y;
			_E_x = E_x;
			throw runtime_error("Unimplemented efield   set_global_params.");
		}

		void set_face_params_facewise(const vector<int> fids, const vector<params_type>& params, bool verbose) override {
			int max_fid = 0;
			for (const auto& fid : fids) {
				if (fid > max_fid) max_fid = fid;
			}
			
			if (max_fid >= _force_enabled_mask_by_cell_index.size()) {
				_force_enabled_mask_by_face_index.resize(max_fid + 1, false);
			}
			
			for (size_t i = 0; i < fids.size(); ++i) {
				int fid = fids.at(i);
				const params_type& fparam = params.at(i);
				
				_force_enabled_mask_by_face_index.at(fid) = true;
			}
		}
    
        bool enabled_for_vertexidx(int fid, bool verbose) {
            if (fid >= _force_enabled_mask_by_face_index.size()) return false;
            
            bool is_enabled = _force_enabled_mask_by_face_index[vid];
            
            if (verbose && is_enabled) {
                cout << "            Confirmed that efield force is on for fid=" << fid << endl;
            }
            
            return is_enabled;
        }
    
    private:
        // vector<Vec> _force_by_vidx;
		float _E_x;
		float _E_y;
		
        vector<bool> _force_enabled_mask_by_face_index;
	};

}

#endif

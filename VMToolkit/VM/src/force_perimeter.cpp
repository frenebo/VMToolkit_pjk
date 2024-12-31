/*!
 * \file force_perimeter.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 30-Nov-2023
 * \brief ForcePerimeter class 
*/ 

#include "force_perimeter.hpp"

using std::cout;
using std::endl;

namespace VMTutorial
{
  
	void ForcePerimeter::compute_all_vertex_forces(vector<Vec>& res, bool verbose)
	{
    
		if (verbose) {
			cout << "   ForcePerimeter::compute_all_vertex_forces - starting" << endl;
    }
		
		size_t n_vertices = _sys.cmesh().cvertices().size();	
		res.resize(n_vertices, Vec(0.0,0.0));
			
		for (auto& vertex : _sys.cmesh().cvertices())
		{
      Vec v_force(0.0,0.0);
			for (auto& he : vertex.circulator()) {
				if (verbose) {
					cout << "    computing force on vertex " << vertex.id << " by halfedge " << he.idx() << endl;
				}
        
				v_force += _compute_he_force(vertex, he, verbose);
			}
      res.at(vertex.id) = v_force;
		}
	}
  
  Vec ForcePerimeter::_compute_he_force(const Vertex& v, const HalfEdge& he, bool verbose)
  {
    if (verbose)
    {
      cout << "        ForcePerimeter::compute - computing force for vertex " << v.id << ", halfedge idx " << he.idx()  << endl;
    }
    
    Vec l = he.to()->data().r - v.data().r;                    // vector along the junction pointing away from the vertex
    const Face& f   = *(he.face());         // cell to the right of the half edge
    const Face& fp = *(he.pair()->face()); // pair cell (opposite side of the same junction)
    
    bool enabled_for_f = enabled_for_faceidx(f.id, verbose);
    bool enabled_for_fp = enabled_for_faceidx(fp.id, verbose);
    
    if ((!enabled_for_f) && (!enabled_for_fp)) {
      return Vec(0.0,0.0);
    }
    
    double P1 = _sys.cmesh().perim(f);
    double P2 = _sys.cmesh().perim(fp);
    
    double gamma_1 = 0.0;
    double lambda_1 = 0.0;
    if (enabled_for_f)  {
      gamma_1 = (f.outer)    ? 0.0 : _gamma.at(f.id);
      lambda_1 = (f.outer)   ? 0.0 : _lambda.at(f.id);
    }
    
    double gamma_2 = 0.0;
    double lambda_2 = 0.0;
    if (enabled_for_fp) {
      gamma_2 = (fp.outer)   ? 0.0 : _gamma.at(fp.id);
      lambda_2 = (fp.outer)  ? 0.0 : _lambda.at(fp.id);
    }

    double lambda = lambda_1 + lambda_2;
    double fedges = gamma_1*P1 + gamma_2*P2 - lambda;
    Vec  perim_force = fedges*l.unit();
		if (verbose)
		{
			cout << "            perim force is (" << perim_force.x << ", " << perim_force.y << ")" << endl;
		}

    return perim_force;
  }
  
  void ForcePerimeter::set_face_params_facewise(const vector<int>& fids, const vector<params_type>& params, bool verbose)
  {
			// Expand the params if necessary
			int max_fid = 0;
			for (const auto& fid : fids) {
			if (fid > max_fid) max_fid = fid;
			}
			
			if (max_fid >= _force_enabled_mask_by_cell_index.size()) {
				_force_enabled_mask_by_cell_index.resize(max_fid + 1, false);
				_gamma.resize(max_fid + 1, 0.0);
				_lambda.resize(max_fid + 1, 0.0);
			}
			
			// Set the param values in the vectors
			for (size_t i = 0; i < fids.size(); ++i) {
				int fid = fids[i];
				const params_type& fparam = params.at(i);
				
				_force_enabled_mask_by_cell_index.at(fid) = true;
				_gamma.at(fid) = fparam.at("gamma");
				_lambda.at(fid) = fparam.at("lambda");
			}
		}

  bool ForcePerimeter::enabled_for_faceidx(int fid, bool verbose) {
    if (fid >= _force_enabled_mask_by_cell_index.size()) return false;
    
    bool is_enabled = _force_enabled_mask_by_cell_index[fid];
    
    if (verbose)
    {
      if (is_enabled) {
        cout << "            Confirmed that perimeter force is on for fid=" << fid << endl;
      }
    }
    
    return is_enabled;
  }
}

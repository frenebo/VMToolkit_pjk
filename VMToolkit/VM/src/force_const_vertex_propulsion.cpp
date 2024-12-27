/*!
 * \file force_const_vertex_propulsion.cpp
 * \author Paul Kreymborg, frenebo@gmail.com
 * \date 17-Dec-2024
 * \brief ForceConstVertexPropulsion class 
*/ 

#include "force_const_vertex_propulsion.hpp"

namespace VMTutorial
{
  void ForceConstVertexPropulsion::set_vertex_params_vertexwise(const vector<int>& vids, const vector<params_type>& params, bool verbose)
  {
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
  
  bool ForceConstVertexPropulsion::enabled_for_vertexidx(int vid, bool verbose) {
    if (vid >= _force_enabled_mask_by_vertex_index.size()) return false;
    
    bool is_enabled = _force_enabled_mask_by_vertex_index[vid];
    
    if (verbose && is_enabled) {
        cout << "            Confirmed that constant propulsion force is on for vid=" << vid << endl;
    }
    
    return is_enabled;
  }
  
  Vec  ForceConstVertexPropulsion::compute_he_force(const Vertex<Property>& v, const HalfEdge<Property>& he, bool verbose)
  {
    if (verbose)
    {
      cout << "        ForceConstVertexPropulsion::compute - computing force for vertex " << v.id << ", halfedge idx " << he.idx()  << endl;
    }

    if (enabled_for_vertexidx(v.id, verbose))
    {
        // This will get called  for each half edge, so we divide by the coordination (aka num of connected vertices)
        Vec fprop_tot = _force_by_vidx.at(v.id);
        
        Vec fprop_divided = Vec(
            fprop_tot.x / v.coordination,
            fprop_tot.y / v.coordination
        );
        if (verbose) {
            cout << "            Force is " << fprop_divided.x <<", " << fprop_divided.y << endl;
        }
        
        return fprop_divided;
    }
    else
    {
        return Vec(0.0,0.0);
    }
  }
}

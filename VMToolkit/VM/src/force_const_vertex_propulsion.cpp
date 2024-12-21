/*!
 * \file force_const_vertex_propulsion.cpp
 * \author Paul Kreymborg, frenebo@gmail.com
 * \date 17-Dec-2024
 * \brief ForceConstVertexPropulsion class 
*/ 

#include "force_const_vertex_propulsion.hpp"

namespace VMTutorial
{
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

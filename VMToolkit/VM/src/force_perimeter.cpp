/*!
 * \file force_perimeter.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 30-Nov-2023
 * \brief ForcePerimeter class 
*/ 

#include "force_perimeter.hpp"

namespace VMTutorial
{
  Vec ForcePerimeter::compute(const Vertex<Property>& v, const HalfEdge<Property>& he, bool verbose)
  {
    if (verbose)
    {
      cout << "ForcePerimeter::compute - computing force for vertex " << v.id << ", halfedge idx " << he.idx()  << endl;
    }
    
    Vec l = he.to()->r - v.r;                    // vector along the junction pointing away from the vertex
    const Face<Property>& f   = *(he.face());         // cell to the right of the half edge
    const Face<Property>& fp = *(he.pair()->face()); // pair cell (opposite side of the same junction)
    double P1 = _sys.mesh().perim(f);
    double P2 = _sys.mesh().perim(fp);
    
    double gamma_1 = 0.0;
    double lambda_1 = 0.0;
    if (enabled_for_faceidx(f.id, verbose))  {
      gamma_1 = (f.outer)    ? 0.0 : _gamma.at(f.id);
      lambda_1 = (f.outer)   ? 0.0 : _lambda.at(f.id);
    }
    
    double gamma_2 = 0.0;
    double lambda_2 = 0.0;
    if (enabled_for_faceidx(fp.id, verbose)) {
      gamma_2 = (fp.outer)   ? 0.0 : _gamma.at(fp.id);
      lambda_2 = (fp.outer)  ? 0.0 : _lambda.at(fp.id);
    }

    double lambda = lambda_1 + lambda_2;
    double fedges = gamma_1*P1 + gamma_2*P2 - lambda;
    Vec  perim_force = fedges*l.unit();
		if (verbose)
		{
			cout << "    perim force is (" << perim_force.x << ", " << perim_force.y << ")" << endl;
		}

    return perim_force;
  }

  double ForcePerimeter::tension(const HalfEdge<Property>& he)
  {
    const Face<Property>& f   = *(he.face());         // cell to the right of the half edge
    const Face<Property>& fp  = *(he.pair()->face()); // pair cell (opposite side of the same junction)
    
    // double lambda_1, lambda_2;
    
    double gamma_1 = 0.0;
    double lambda_1 = 0.0;
    if (enabled_for_faceidx(f.id, false)) {
      gamma_1 = (f.outer)    ? 0.0 : _gamma.at(f.id);
      lambda_1 = (f.outer)   ? 0.0 : _lambda.at(f.id);
    }
    
    double gamma_2 = 0.0;
    double lambda_2 = 0.0;
    if (enabled_for_faceidx(fp.id, false)) {
      gamma_2 = (fp.outer)   ? 0.0 : _gamma.at(fp.id);
      lambda_2 = (fp.outer)  ? 0.0 : _lambda.at(fp.id);
    }

    double lambda = lambda_1 + lambda_2;
  
    return gamma_1*_sys.mesh().perim(f) + gamma_2*_sys.mesh().perim(fp) - lambda;
  }

  double ForcePerimeter::energy(const Face<Property>& f)
  {
    if (f.outer || f.erased)
      return 0.0;

    double P = _sys.mesh().perim(f);
    
    
    double gamma = 0.0;
    double lambda = 0.0;
    
    if (enabled_for_faceidx(f.id, false)) {
      gamma = _gamma.at(f.id);
      lambda = _lambda.at(f.id);
    }
    
    
    double P0;
    
    P0 = lambda/gamma;
    
    double dP = P - P0;

    return 0.5*gamma*dP*dP;
    
  }

}

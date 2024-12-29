/*!
 * \file integrator.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com 
 * \date 30-Nov-2023
 * \brief IntegratorEuler class 
*/

#include "integrator_euler.hpp"

namespace VMTutorial
{
  void IntegratorEuler::step(bool verbose)
  {
    if (verbose) {
      cout << "IntegratorEuler::step" << endl;
    }
    
    double mu = 1.0 / _gamma;    // mobility 
    double B = sqrt(2.0*mu*_T);
    double sqrt_dt = sqrt(_dt);
    
    // Compute force on each vertex
    if (verbose) {
      cout << "Total number of vertices in system - " << _sys.mesh().vertices().size() << endl;
    }
    
    std::vector<Vec> vertex_forces = _force_compute.compute_all_vertex_forces(verbose);
    
    if (verbose) {
      cout << "Done computing vertex forces" << endl;
    }
    
    // This is actual integrator 
    for (size_t vid = 0; vid < _sys.mesh().vertices().size(); vid++)
    {
      Vertex<Property>& v = _sys.mesh().vertices()[vid];
      
      if (!v.erased)
      {
        const Vec& f = vertex_forces.at(vid);

        if (verbose) {
          cout << "     Adding to the r vector" << endl;
        }
        Vec rold = v.data().r;
        v.data().r += _dt*mu*f;
        
        if (verbose) {
          cout << "     added _dt*mu*f" << endl;
        }
        
        // if (_T > 0.0)
        // {
        //   Vec ffr(B*_rng.gauss_rng(), B*_rng.gauss_rng());  // random noise contribution to force
        //   Vec fr = ffr;
          
        //   v.data().r += sqrt_dt*fr;  // update vertex position due to noise
        // }
        v.data().vel = (1.0 / _dt) * (v.data().r - rold);  
      }
    }
    if (verbose) {
      cout << "Finished second round thingy" << endl;
    }
  }

}
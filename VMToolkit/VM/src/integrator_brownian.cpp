/*!
 * \file integrator.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com 
 * \date 30-Nov-2023
 * \brief IntegratorBrownian class 
*/

#include "integrator_brownian.hpp"

namespace VMTutorial
{
  void IntegratorBrownian::step(bool verbose)
  {
    // bool verbose = true;
    if (verbose) {
      cout << "IntegratorBrownian::step" << endl;
    }
    // cout << "Doing"
    double mu = 1.0 / _gamma;    // mobility 
    double B = sqrt(2.0*mu*_T);
    double sqrt_dt = sqrt(_dt);
    // Compute force on each vertex
    size_t logvi = 0;
    if (verbose) {
      cout << "Total number of vertices in system - " << _sys.mesh().vertices().size() << endl;
    }
    _force_compute.compute_and_set_all_vertex_forces(verbose);
    // for (auto& v : _sys.mesh().vertices())
    // {
    //   if (!v.erased)
    //   {
    //     if (verbose) {
    //       cout << "  first round Computing force on vertex " << logvi << endl;
    //     }
    //     logvi++;
    //     v.data().force = _force_compute.compute_v_force(v, verbose);
    //   }
    // }
    if (verbose) {
      cout << "Done computing vertex forces" << endl;
    }
    
    // This is actual integrator 
    logvi = 0;
    for (auto& v : _sys.mesh().vertices())
    {
      // if (verbose) {
      //   cout << "   second round " << logvi++ << endl;
      // }
      if (!v.erased)
      {
        // if (verbose) {
        //   cout << "     adding force to constant force" << endl;
        // }
        // if (v)
        // add external force 
        // if (verbose) {
        //   // cout << "vert type is " << v.data().vert_type << endl;
        //   // cout << "_constant_force size: " << _constant_force.size() << endl;
        // }
        
        // @TODO check that the _constant_force has been properly set up at this point.
        Vec f = v.data().force;
        
        // // Add external force if present for this vertex id...
        // map<int,Vec>::iterator extf_it = _const_ext_forces_by_vid.find(v.id);
        // // Bar b3;
        // if(extf_it != _const_ext_forces_by_vid.end())
        // {
        //   //element found;
        //   Vec ext_force_on_vtx = extf_it->second;
        //   f += ext_force_on_vtx;
        // }
        
        // apply constraint
        if (_constraint_enabled) {
          if (verbose) {
            cout << "     Applying constraint " << endl;
          }
          f = _constrainer->apply_vertex(v, f);
        }

        if (verbose) {
          cout << "     Adding to the r vector" << endl;
        }
        Vec rold = v.r;
        v.r += _dt*mu*f;  // deterministic part of the integrator step
        
        if (verbose) {
          cout << "     added _dt*mu*f" << endl;
        }
        
        if (_T > 0.0)
        {
          Vec ffr(B*_rng.gauss_rng(), B*_rng.gauss_rng());  // random noise contribution to force
          Vec fr = ffr;
          if (_constraint_enabled) 
            fr = _constrainer->apply_vector(v, ffr);
          v.r += sqrt_dt*fr;  // update vertex position due to noise
        }
        v.data().vel = (1.0 / _dt) * (v.r - rold);  
      }
    }
    if (verbose) {
      cout << "Finished second round thingy" << endl;
    }

    // Update direction of cell director using simple Brownian dynamics
    if (_update_n)
    {
      double stoch_coeff = sqrt(_Dr*_dt);
      for (auto& f : _sys.mesh().faces())
      {
        if (!(f.erased || f.outer))
        {
          double dtheta = stoch_coeff*_rng.gauss_rng();
          Vec n = f.data().n;
          // Sine and cosine of the rotation angle
          double c = cos(dtheta), s = sin(dtheta);
          // Rotation matrix around z axis
          double Rxx = c, Rxy = -s;
          double Ryx = s, Ryy = c;
          // Apply rotation matrix
          double nx = Rxx*n.x + Rxy*n.y;
          double ny = Ryx*n.x + Ryy*n.y;
          double len = sqrt(nx*nx + ny*ny);
          f.data().n = Vec(nx/len, ny/len);  // Rotation does not change length of a vector, but numerical errors can accumulate, so we normalize it
        }
      }
    }

  }

}
/*!
 * \file integrator_adaptive_runge_kutta.cpp
 * \author Paul Kreymborg, frenebo@gmail.com
 * \date 31-Dec-2024
 * \brief IntegratorAdaptiveRungeKutta class 
*/

#include "integrator_adaptive_runge_kutta.hpp"


namespace VMTutorial
{
  vector<Vec> IntegratorAdaptiveRungeKutta::_instantaneous_velocities(bool verbose) const
  {
    double mu = 1.0 / _gamma;    // mobility 
    if (verbose) {
      cout << "  calculating instantaneous velocities with mu=1.0/_gamma=" << mu << endl;
    }
    
    vector<Vec> vertex_vels;
    _force_compute.compute_all_vertex_forces(vertex_vels, verbose);
    
    // v = F/gamma
    for (auto& vel : vertex_vels)
    {
      vel.x *= mu;
      vel.y *= mu;
    }
    
    return vertex_vels;
  }
  
  void IntegratorAdaptiveRungeKutta::step(bool verbose)
  {
    if (_gamma_param <= 0) {
      throw runtime_error(" IntegratorAdaptiveRungeKutta::step - 'gamma' has not been set to valid value");
    }
    if (_error_allowed_param <= 0) {
      throw runtime_error(" IntegratorAdaptiveRungeKutta::step - 'error_alloweds' has not been set to valid value");
    }
    
    throw runtime_error("Not implemented - IntegratorAdaptiveRungeKutta::step ");
  }
}
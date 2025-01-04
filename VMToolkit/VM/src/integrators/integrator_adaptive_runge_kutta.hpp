/*!
 * \file integrator_adaptive_runge_kutta.hpp
 * \author Paul Kreymborg, frenebo@gmail.com
 * \date 31-Dec-2024
 * \brief IntegratorAdaptiveRungeKutta class 
*/

#ifndef __INTEGRATOR_ADAPTIVE_RUNGE_KUTTA_HPP__
#define __INTEGRATOR_ADAPTIVE_RUNGE_KUTTA_HPP__


#include "integrator.hpp"

#include <map>
#include <stdexcept>

using std::map;
using std::runtime_error;


namespace VMTutorial
{

  class IntegratorAdaptiveRungeKutta : public Integrator 
  {
    /*
      Implementation of RK4
    */
    public:
      
      IntegratorAdaptiveRungeKutta(
        System& sys,
        ForceCompute& fc,
        int seed
      ) : Integrator{sys, fc, seed},
          _gamma{-1.0},
          _error_allowed{-1.0}
      {
      }
      
      void step(bool verbose) override;
      
      void set_params(const params_type& params) override 
      { 
        for (auto& p : params)
        {
          if (p.first == "gamma")
          {
            _gamma = p.second;
          } else if (p.first == "error_allowed") {
            _error_allowed = p.second;
          } else {
            throw runtime_error("Unknown parameter name - " + p.first);
          }
        }
      }
      
      
    private:
      vector<Vec> _instantaneous_velocities(bool verbose) const;
      
      double _gamma_param;
      double _error_allowed_param;
  };
}


#endif

/*!
 * \file integrator_runge_kutta.hpp
 * \author Paul Kreymborg, frenebo@gmail.com
 * \date 27-Dec-2024
 * \brief IntegratorRungeKutta class 
*/

#ifndef __INTEGRATOR_RUNGE_KUTTA_HPP__
#define __INTEGRATOR_RUNGE_KUTTA_HPP__


#include "integrator.hpp"

#include <map>
#include <stdexcept>

using std::map;
using std::runtime_error;

namespace VMSim
{

  class IntegratorRungeKutta : public Integrator 
  {
    /*
      Implementation of RK4
    */
    public:
      
      IntegratorRungeKutta(
        System& sys,
        ForceCompute& fc,
        Topology& top
      ) : Integrator{
          sys,
          fc,
          top,
          Integrator::IntegrationType::MANUAL_TIMESTEP
        },
        _gamma{-1.0},
        _dt{-1.0}
      {
      }
      
      void timestep_manual(bool verbose) override;
      
      void set_params(const params_type& params) override 
      { 
        for (auto& p : params)
        {
          if (p.first == "gamma")
          {
            _gamma = p.second;
          } else if (p.first == "dt") {
            _dt = p.second;
          } else {
            throw runtime_error("Unknown parameter name - " + p.first);
          }
        }
      }
      
      
    private:
      vector<Vec> _instantaneous_velocities(bool verbose) const;
      
      double _dt;
      double _gamma;             // friction 
  };
}


#endif

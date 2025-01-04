/*!
 * \file integrator_euler.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 30-Nov-2023
 * \brief IntegratorEuler class 
*/

#ifndef __INTEGRATOR_EULER_HPP__
#define __INTEGRATOR_EULER_HPP__


#include "integrator.hpp"

#include <map>
#include <string>

using std::map;
using std::runtime_error;

namespace VMTutorial
{

  class IntegratorEuler : public Integrator 
  {

    public:

      IntegratorEuler(System& sys, ForceCompute& fc, int seed) : Integrator{sys, fc, seed},
                                                                    _T{0.0},
                                                                    _gamma{1.0},
                                                                    _Dr{0.0},
                                                                    _dt{1.0}
                                                                    
      {
      }

      void step(bool verbose) override;
      void set_params(const params_type& params) override 
      { 
        for (auto& p : params)
        {
          if (p.first == "T") {
            _T = p.second;
          } else if (p.first == "gamma") {
            _gamma = p.second;
          } else if (p.first == "Dr") {
            _Dr = p.second;
          }  else if (p.first == "dt") {
            _dt = p.second;
          } else {
            throw runtime_error("Unknown parameter name - " + p.first);
          }
        }
      }
      
    private:

      double _T;                 // temperature
      double _dt;
      double _gamma;             // friction 
      double _Dr;                // rotational diffusion constant
      // bool _update_n;            // If true, update direction of the cell self-propulsion direction
      
    };

}

#endif


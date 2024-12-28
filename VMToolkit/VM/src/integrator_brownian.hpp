/*!
 * \file integrator_brownian.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 30-Nov-2023
 * \brief IntegratorBrownian class 
*/

#ifndef __INTEGRATOR_BROWNIAN_HPP__
#define __INTEGRATOR_BROWNIAN_HPP__

// #include "constrainer.hpp"
// #include "constraint_none.hpp"
// #include "constraint_fixed.hpp"


#include "integrator.hpp"

#include <chrono>
#include <utility>
#include <map>
#include <string> 
#include <memory>

using namespace std::chrono;
using std::map;
using std::make_unique;

namespace VMTutorial
{

  // using ConstrainerType = unique_ptr<Constrainer>;

  class IntegratorBrownian : public Integrator 
  {

    public:

      IntegratorBrownian(System& sys, ForceCompute& fc, int seed) : Integrator{sys, fc, seed},
                                                                    _T{0.0},
                                                                    _gamma{1.0},
                                                                    _Dr{0.0},
                                                                    _update_n{false}
                                                                    
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
          } else {
            throw runtime_error("Unknown parameter name - " + p.first);
          }
        }
      }
      
      
      void set_flag(const string& flag) override 
      {  
        if (flag == "update_n")
          _update_n = true;
        else
          throw runtime_error("Brownian integrator: Unknown flag : " + flag + ".");
      }

      
    private:

      double _T;                 // temperature
      double _gamma;             // friction 
      double _Dr;                // rotational diffusion constant
      bool _update_n;            // If true, update direction of the cell self-propulsion direction
      
    };

}

#endif


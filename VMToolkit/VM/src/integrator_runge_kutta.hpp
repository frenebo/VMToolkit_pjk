/*!
 * \file integrator_runge_kutta.hpp
 * \author Paul Kreymborg, frenebo@gmail.com
 * \date 27-Dec-2024
 * \brief IntegratorRungeKutta class 
*/

#ifndef __INTEGRATOR_RUNGE_KUTTA_HPP__
#define __INTEGRATOR_RUNGE_KUTTA_HPP__


#include "constrainer.hpp"
#include "constraint_none.hpp"
#include "constraint_fixed.hpp"
#include "integrator.hpp"

using namespace std::chrono;
using std::map;
using std::make_unique;

namespace VMTutorial


namespace VMTutorial
{

  using ConstrainerType = unique_ptr<Constrainer>;


  class IntegratorRungeKutta : public Integrator 
  {
    public:
      
      IntegratorRungeKutta(
        System& sys,
        ForceCompute& fc,
        int seed,
      ) : Integrator{sys, fc, seed},
          _gamma{1.0}
      {
        _constrainer = make_unique<Constrainer>();
        _constrainer->add<ConstraintNone>("none");
        _constrainer->add<ConstraintFixed>("fixed");
      }
      
      void step(bool verbose) override;
      
      void set_params(const params_type& params) override 
      { 
        for (auto& p : params)
        {
          if (p.first == "gamma")
          {
            _gamma = p.second;
          } else {
            throw runtime_error("Unknown parameter name - " + p.first);
          }
        }
      }
      
      
      void set_string_params(const string_params_type& params) override {
        throw runtime_error("IntegratorRungeKutta::set_string_params - unimplemented");
      }
      
      
      void set_flag(const string& flag) override 
      {  
          throw runtime_error("IntegratorRungeKutta::set_flag - Unknown flag : " + flag + ".");
      }
      
    private:
      ConstrainerType _constrainer; // Apply various constraints
      double _gamma;             // friction 
  }
  {


#endif

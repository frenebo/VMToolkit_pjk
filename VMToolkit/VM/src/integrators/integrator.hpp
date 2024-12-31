/*!
 * \file integrator.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 30-Nov-2023
 * \brief Integrator class 
*/

#ifndef __INTEGRATOR_HPP__
#define __INTEGRATOR_HPP__

#include "../rng.hpp"
#include "../system.hpp"

#include "../force_compute.hpp"

#include <stdexcept>
#include <chrono>
#include <map>

using std::map;
using std::runtime_error;

namespace VMTutorial
{

  class Integrator 
  {

    public:
    
      Integrator(System& sys, ForceCompute& fc, int seed) : _sys{sys}, 
                                                            _force_compute{fc}, 
                                                            _rng{RNG((seed >= 0) ? seed : std::chrono::system_clock::now().time_since_epoch().count())},
                                                            _dt{0.01}
      { 
        
      }
      virtual ~Integrator() { }
      
      virtual void step(bool verbose=false)
      {
        throw runtime_error("Child has not implemented Integrator::step");
      }
      
      virtual void set_params(const params_type&)
      {
        throw runtime_error("Child has not implemented Integrator::set_params");
      }
      
      void set_dt(double dt) { _dt = dt; }
      void rng_set(const RNGState& state) { _rng.set(state); }
      RNGState get_rng_state() { return _rng.get_state(); }

    protected:

      RNG _rng;
      System& _sys;              // system
      ForceCompute&   _force_compute;
      double _dt; // time step
  };

}

#endif


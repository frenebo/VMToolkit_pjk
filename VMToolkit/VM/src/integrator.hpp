/*!
 * \file integrator.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 30-Nov-2023
 * \brief Integrator class 
*/

#ifndef __INTEGRATOR_HPP__
#define __INTEGRATOR_HPP__

#include "rng.hpp"
#include "system.hpp"
#include "force.hpp"

#include "constrainer.hpp"
#include "constraint_none.hpp"
#include "constraint_fixed.hpp"

#include "force_compute.hpp"

#include <chrono>
#include <utility>
#include <map>

using namespace std::chrono;
using std::map;

namespace VMTutorial
{

  class Integrator 
  {

    public:
    
      Integrator(System& sys, ForceCompute& fc, int seed) : _sys{sys}, 
                                                            _force_compute{fc}, 
                                                            _rng{RNG((seed >= 0) ? seed : system_clock::now().time_since_epoch().count())},
                                                            _enabled{true},
                                                            _constraint_enabled{true},
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
      // virtual void set_params_for_vertices_elementwise(const vector<int>& vids, const vector<params_type>& params) = 0;
      // virtual void set_params_for_faces_elementwise(const vector,int>& fids, const vector<params_type>& params) = 0;
      
      virtual void set_string_params(const string_params_type& params)
      {
        throw runtime_error("Child has not implemented Integrator::set_string_params");
      }
      // // virtual void set_external_force(const string&, const Vec&) = 0;
      // virtual void set_external_forces_by_vertex(const vector<int>& vids, const vector<Vec>& forces)
      // {
      //   throw runtime_error("Child has not implemented Integrator::set_external_forces_by_vertex");
      // }
      virtual void set_flag(const string&)
      {
        throw runtime_error("Child has not implemented Integrator::set_flag");
      }

      void set_dt(double dt) { _dt = dt; }
      void enable() { _enabled = true; }
      void disable() { _enabled = false; }
      bool is_enabled() { return _enabled; }
      void enable_constraint(bool enable) { _constraint_enabled = enable;  }
      void rng_set(const RNGState& state) { _rng.set(state); }
      RNGState get_rng_state() { return _rng.get_state(); }

    protected:

      RNG _rng;
      System& _sys;              // system
      ForceCompute&   _force_compute; 
      bool _enabled;
      bool _constraint_enabled;
      double _dt; // time step
  };

}

#endif


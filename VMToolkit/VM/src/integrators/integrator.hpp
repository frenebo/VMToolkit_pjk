/*!
 * \file integrator.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 30-Nov-2023
 * \brief Integrator class 
*/

#ifndef __INTEGRATOR_HPP__
#define __INTEGRATOR_HPP__

// #include "../rng.hpp"
#include "../system.hpp"
#include "../force_compute.hpp"
#include "../topology.hpp"

#include <stdexcept>
#include <chrono>
#include <map>


namespace VMSim
{
  using std::map;
  using std::runtime_error;

  class Integrator 
  {
    public:
      enum class IntegrationType { MANUAL_TIMESTEP, ADAPTIVE_TIMESTEP };

    
      Integrator(
        System& sys,
        ForceCompute& fc,
        Topology& top,
        IntegrationType integ_type
      ) : _sys{sys}, 
          _force_compute{fc},
          _topology{top},
          _integ_type{integ_type}
      { 
        
      }
      virtual ~Integrator() { }
      
      virtual void timestep_manual(bool verbose)
      {
        if (_integ_type == IntegrationType::ADAPTIVE_TIMESTEP) {
          throw runtime_error("Integrator::timestep_manual - this is an ADAPTIVE_TIMESTEP integrator, should not be run with timestep_manual");
        }
        
        throw runtime_error("Child has not implemented Integrator::timestep_manual");
      }
      
      virtual void adaptive_run_for_time(double t_run, bool show_progress_bar, bool topological_change, bool verbose)
      {
        if (_integ_type == IntegrationType::MANUAL_TIMESTEP) {
          throw runtime_error("Integrator::adaptive_run_for_time - this is a MANUAL_TIMESTEP integrator, should not be run with adaptive_run_for_time");
        }
        
        throw runtime_error("Child has not implemented Integrator::adaptive_run_for_time");
      }
      
      virtual void set_params(const params_type&)
      {
        throw runtime_error("Child has not implemented Integrator::set_params");
      }
      
      virtual void set_flag(const string& flagname, bool flag_val)
      {
        throw runtime_error("Integrator::set_flag - child has not implemented set_flag");
      }
      
      IntegrationType get_integration_type()
      {
        return _integ_type;
      }
      
    protected:
      System& _sys;              // system
      ForceCompute&   _force_compute;
      Topology& _topology;
      IntegrationType _integ_type;
  };

}

#endif


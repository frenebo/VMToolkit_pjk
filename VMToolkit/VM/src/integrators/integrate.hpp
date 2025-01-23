/*!
 * \file integrate.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 30-Nov-2023
 * \brief Integrate class 
 */
 
#ifndef __INTEGRATE_COMPUTE_HPP__
#define __INTEGRATE_COMPUTE_HPP__

#include <vector>
#include <string>

#include "../system.hpp"
#include "../force_compute.hpp"
#include "integrator.hpp"


namespace VMSim
{
  using std::runtime_error;
  using std::vector;
  using std::string;
  using std::map;
  using std::endl;
  using std::cout;

  class Integrate : public ClassFactory<Integrator>
  {
    public:

      Integrate(System& sys, ForceCompute& fc, Topology& top) : _sys{sys},
                                                _force_compute{fc},
                                                _topology{top}
      { 

      }
      
      ~Integrate() = default; 

      Integrate(const ForceCompute&) = delete;

      void timestep_manual(bool verbose)
      {
        if (verbose) {
          cout << "In Integrate::timestep_manual" << endl;
        }
        
        for (const string& i : this->_integ_order) {
          if (verbose) {
            cout << "Going to apply a step for " << i << endl;
          }
          
          this->factory_map[i]->timestep_manual(verbose);
        }
      }
      
      void run_time_adaptive(double runtime_tot, bool show_progress_bar, bool topological_change, bool verbose)
      {
        if (verbose) {
          cout << "In Integrate::run_time_adaptive - running with runtime_tot=" << runtime_tot  << endl;
        }
      
      
        for (const string& i : this->_integ_order) {
          if (verbose) {
            cout << "Going to apply adaptive runtime for "<< i << endl;
          }
          
          this->factory_map[i]->adaptive_run_for_time(runtime_tot, show_progress_bar, topological_change, verbose);
        }
        
      }
      
      void set_flag(const string& iname, const string& flagname, bool flag_val)
      {
        if (this->factory_map.find(iname) != this->factory_map.end()) {
          this->factory_map[iname]->set_flag(flagname, flag_val);
        } else {
          throw runtime_error("Integrate::set_flag - integrator with name could not be found: " + iname);
        }
      }

      void set_params(const string& iname, const params_type& params)
      {
        if (this->factory_map.find(iname) != this->factory_map.end())
          this->factory_map[iname]->set_params(params);
        else
          throw runtime_error("set_params: Integrator type " + iname + " is not used in this simulation.");
      }
      
      void add_integrator(const string& iname);

    private: 

      System& _sys;
      ForceCompute& _force_compute;
      Topology& _topology;
      // double _dt;
      vector<string> _integ_order;  // Keeps track in which order integrators were added
      map<string,bool> _integrators_enabled;
  };

  void export_Integrate(py::module&);

}
#endif


  

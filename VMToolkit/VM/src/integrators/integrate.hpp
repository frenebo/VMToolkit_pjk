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

using std::runtime_error;
using std::vector;
using std::string;
using std::map;


namespace VMTutorial
{
  class Integrate : public ClassFactory<Integrator>
  {
    public:

      Integrate(System& sys, ForceCompute& fc) : _sys{sys},
                                                _force_compute{fc}
      { 

      }
      
      ~Integrate() = default; 

      Integrate(const ForceCompute&) = delete;

      void timestep_manual(bool verbose=false)
      {
        if (verbose) {
          std::cout << "In Integrate::timestep_manual" << std::endl;
        }
        // int n_integrators = 0;
        for (const string& i : this->_integ_order) {
          // if (this->_integrators_enabled.at(i)) {
          if (verbose) {
            std::cout << "Going to apply a step for " << i << std::endl;
          }
          this->factory_map[i]->timestep_manual(verbose);
          // }
          
          // n_integrators++;
        }
      }

      void set_params(const string& iname, const params_type& params)
      {
        if (this->factory_map.find(iname) != this->factory_map.end())
          this->factory_map[iname]->set_params(params);
        else
          throw runtime_error("set_params: Integrator type " + iname + " is not used in this simulation.");
      }

      // void enable(const string& iname)
      // {
      //   if (this->factory_map.find(iname) != this->factory_map.end()) {
      //     this->_integrators_enabled[iname] = true;
      //   } else {
      //     throw runtime_error("enable: Integrator type " + iname + " is not used in this simulation."); 
      //   }
      // }

      // void disable(const string& iname)
      // {
      //   if (this->factory_map.find(iname) != this->factory_map.end())
      //     this->_integrators_enabled[iname] = false;
      //   else
      //     throw runtime_error("disable: Integrator type " + iname + " is not used in this simulation.");
      // }
      
      void add_integrator(const string& iname);

    private: 

      System& _sys;
      ForceCompute& _force_compute;
      // double _dt;
      vector<string> _integ_order;  // Keeps track in which order integrators were added
      map<string,bool> _integrators_enabled;
  };

  void export_Integrate(py::module&);

}
#endif


  

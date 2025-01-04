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

      Integrate(System& sys, ForceCompute& fc, int seed) : _sys{sys},
                                                           _force_compute{fc},
                                                          //  _dt{0.01}, 
                                                           _seed{seed}
                                                           { 
          
                                                           } 
      ~Integrate() = default; 

      Integrate(const ForceCompute&) = delete;

      void apply(bool verbose=false)
      {
        if (verbose) {
          std::cout << "In Integrate::apply" << std::endl;
        }
        for (const string& i : this->_integ_order) {
          if (this->_integrators_enabled.at(i)) {
            if (verbose) {
              std::cout << "Going to apply a step for " << i << std::endl;
            }
            this->factory_map[i]->step(verbose);
          }
        }
        _sys.time_step()++;
        _sys.simulation_time() += _dt;
      }

      void set_params(const string& iname, const params_type& params)
      {
        if (this->factory_map.find(iname) != this->factory_map.end())
          this->factory_map[iname]->set_params(params);
        else
          throw runtime_error("set_params: Integrator type " + iname + " is not used in this simulation.");
      }

      void enable(const string& iname)
      {
        if (this->factory_map.find(iname) != this->factory_map.end()) {
          this->_integrators_enabled[iname] = true;
        } else {
          throw runtime_error("enable: Integrator type " + iname + " is not used in this simulation."); 
        }
      }

      void disable(const string& iname)
      {
        if (this->factory_map.find(iname) != this->factory_map.end())
          this->_integrators_enabled[iname] = false;
        else
          throw runtime_error("disable: Integrator type " + iname + " is not used in this simulation.");
      }

      // void set_dt(double dt) 
      // { 
      //   if (this->factory_map.size() < 1)
      //     throw runtime_error("There are no integrators defined. Time step cannot be changed.");
      //   _dt = dt; 
      //   for (auto& integ : this->factory_map)
      //     integ.second->set_dt(dt);
      // }

      
      void add_integrator(const string& iname);

    private: 

      System& _sys;
      ForceCompute& _force_compute;
      // double _dt;
      int _seed;
      vector<string> _integ_order;  // Keeps track in which order integrators were added
      map<string,bool> _integrators_enabled;
  };

  void export_Integrate(py::module&);

}
#endif


  

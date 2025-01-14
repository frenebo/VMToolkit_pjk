
/*!
 * \file simulation.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 30-Nov-2023
 * \brief Simulation class 
 */ 

#ifndef __SIMULATION_HPP__
#define __SIMULATION_HPP__

#define XSTR(s) STR(s)
#define STR(s) #s

#include <iostream>
#include <string>

#include "system.hpp"
#include "force_compute.hpp"
#include "integrators/integrate.hpp"
#include "topology.hpp"
#include "dump.hpp"
#include "version.hpp"

using std::cout;
using std::string;
using std::endl;

namespace VMSim
{

  class Simulation
  {
  public:
    Simulation(System& sys, Integrate& integ, ForceCompute& f, Topology& t) : _sys{sys}, 
                                                                              _integ{integ}, 
                                                                              _force_compute{f},
                                                                              _topology{t},
                                                                              print_freq{100},
                                                                              bar_width{40},
                                                                              sim_step{0}                                                                              
    { 

    }
    void run_timestep_manual(int steps, bool topological_change, bool old_style, bool verbose);
    void run_time_adaptive(double runtime_tot, bool topological_change, bool verbose);
    const string print_version() {
      return  "branch : "+static_cast<string>(XSTR(GIT_BRANCH))+" commit : "+static_cast<string>(XSTR(GIT_COMMIT_HASH));
     }
    void progress_bar(double, const string&);


    int print_freq;
    int bar_width;
    int sim_step;
    
  private:
    // variables and parameters
    System& _sys;
    Integrate& _integ;
    ForceCompute& _force_compute;
    Topology& _topology;
    
  };

  void export_Simulation(py::module&);

}

#endif

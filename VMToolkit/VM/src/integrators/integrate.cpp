/*!
 * \file integrate.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 22-May-2019
 * \brief Integrate class 
 */ 

#include "integrate.hpp"
#include "integrator_euler.hpp"
#include "integrator_runge_kutta.hpp"

namespace VMTutorial
{
  void Integrate::add_integrator(const string& iname)
  {
    string name = iname; 
    std::transform(name.begin(), name.end(), name.begin(), ::tolower);
    
    if (_integ_order.size() != 0) {
      throw runtime_error("Integrate::add_integrator - Only one integrator can be active at a time");
    }
    
    if (name == "euler") {
      this->add<IntegratorEuler, System&, ForceCompute&, int>(name, _sys, _force_compute, _seed);
    } else if (name == "runge_kutta") {
      this->add<IntegratorRungeKutta, System&, ForceCompute&, int>(name, _sys, _force_compute, _seed);
    } else  {
      throw runtime_error("Unknown integrator type : " + name + ".");
    }
    _integ_order.push_back(name);
    // _integrators_enabled[iname] = true;
  }
  
  void export_Integrate(py::module& m)
  {
    py::class_<Integrate>(m, "Integrate")
      .def(py::init<System&, ForceCompute&, int>())
      .def("set_params", &Integrate::set_params)
      // .def("enable", &Integrate::enable)
      // .def("disable", &Integrate::disable)
      // .def("set_dt", &Integrate::set_dt)
      .def("add", &Integrate::add_integrator)
      ;
  }
}
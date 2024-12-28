/*!
 * \file integrate.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 22-May-2019
 * \brief Integrate class 
 */ 

#include "integrate.hpp"

namespace VMTutorial
{
  void export_Integrate(py::module& m)
  {
    py::class_<Integrate>(m, "Integrate")
      .def(py::init<System&, ForceCompute&, int>())
      .def("set_params", &Integrate::set_params)
      .def("enable", &Integrate::enable)
      .def("disable", &Integrate::disable)
      .def("set_dt", &Integrate::set_dt)
      .def("add", &Integrate::add_integrator)
      ;
  }
}
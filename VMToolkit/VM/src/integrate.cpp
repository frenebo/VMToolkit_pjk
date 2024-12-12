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
      .def("set_string_params", &Integrate::set_string_params)
      .def("set_external_forces_by_vertex", &Integrate::set_external_forces_by_vertex)
      .def("set_flag", &Integrate::set_flag)
      .def("enable", &Integrate::enable)
      .def("disable", &Integrate::disable)
      .def("set_dt", &Integrate::set_dt)
      .def("enable_constraint", &Integrate::enable_constraint)
      .def("add", &Integrate::add_integrator);
  }
}
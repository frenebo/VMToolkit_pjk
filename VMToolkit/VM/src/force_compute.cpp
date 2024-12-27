/*!
 * \file force_compute.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 30-Nov-2023
 * \brief ForceCompute class
 */

#include "force_compute.hpp"

namespace VMTutorial
{
	void export_ForceCompute(py::module &m)
	{
		py::class_<ForceCompute>(m, "ForceCompute")
			.def(py::init<System &>())
			// .def("set_params", &ForceCompute::set_params)
			// .def("set_vec_params", &ForceCompute::set_vec_params)
			.def("set_global_params", &ForceCompute::set_global_params, py::arg("force_id"), py::arg("num_params"), py::arg("str_params"), py::arg("verbose")=false)
			.def("set_face_params_facewise", &ForceCompute::set_face_params_facewise)
			.def("set_vertex_params_vertexwise", &ForceCompute::set_vertex_params_vertexwise)
			// .def("set_flag", &ForceCompute::set_flag)
			.def("add_force", &ForceCompute::add_force, py::arg("force_id"), py::arg("force_type"), py::arg("verbose")=false)
			// .def()

			.def("start_force_compute_timers", &ForceCompute::start_force_compute_timers, py::arg("verbose")=false)
			.def("get_force_compute_timers_millis", &ForceCompute::get_force_compute_timers_millis, py::arg("verbose")=false);
			// .def("compute", &ForceCompute::compute_forces)
			// .def("energy", &ForceCompute::total_energy);
	}
}
/*!
 * \file simulation.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 30-Nov-2023
 * \brief Simulation class 
 */ 

#include "simulation.hpp"

using std::runtime_error;
using std::cout;
using std::endl;
namespace py = pybind11;


namespace VMSim
{
  void Simulation::run_time_adaptive(double runtime_tot, bool topological_change, bool verbose)
  {
    if (verbose) {
      cout << "Simulation::run_time_adaptive - running simulation for runtime_tot=" << runtime_tot << endl;
    }
    
    
    // throw runtime_error("This should integrate T1 transitions into the integrator...");
    
    bool show_progress_bar=true;
    _integ.run_time_adaptive(runtime_tot, show_progress_bar, topological_change, verbose);
  }
  
  void Simulation::run_timestep_manual(int steps, bool topological_change, bool old_style, bool verbose)
  {
    double progress = 0.0;
    if (verbose) {
      cout << "Simulation::run_timestep_manual - Running simulation for " << steps << " steps" << endl;
    }
    
    for (int i = sim_step; i < sim_step + steps; i++)
    {
      if (verbose) { cout << "doing step #" << i << endl; }
      // _sys.set_topology_change(false);
      if (topological_change)
      {
          _topology.T1(verbose);
      }
      
      if (verbose) { cout << "doing integration" << endl; }
      _integ.timestep_manual(verbose);

      if (this->print_freq > 0) 
      {
        if (old_style)
        {
          if (i % this->print_freq == 0)
            cout << "step : " << i << endl;
        }
        else 
        {
          progress_bar(progress, "\r");
          progress += 1.0 / steps;
        }
        
      }
    }
    
    if (this->print_freq > 0 && !old_style) {
      progress_bar(progress, " ");
    }
    
    sim_step += steps;
    
    if (this->print_freq > 0 && !old_style) {
      cout << " --> Completed " << sim_step << " simulation steps " << endl;  
    }
  }
  
  

  void Simulation::progress_bar(double progress, const string& end_of_line)
  {
    cout << "[";
    int pos = static_cast<int>(round(bar_width * progress));
    for (int i = 0; i < bar_width; ++i) 
    {
      if (i < pos) 
        cout << "=";
      else if (i == pos) 
        cout << ">";
      else 
        cout << " ";
    }
    cout << "] "  << static_cast<int>(round(progress * 100.0)) << "%" << end_of_line;
    cout.flush();
  }


  void export_Simulation(py::module& m)
  {
    py::class_<Simulation>(m, "Simulation")
          .def(py::init<System&,Integrate&,ForceCompute&,Topology&>())
          .def_readwrite("print_freq", &Simulation::print_freq)
          .def_readwrite("bar_width", &Simulation::bar_width)
          .def("run_time_adaptive", &Simulation::run_time_adaptive, py::arg("runtime_tot"),  py::arg("topological_change") = true, py::arg("verbose") = false)
          .def("run_timestep_manual", &Simulation::run_timestep_manual, py::arg("steps"), py::arg("topological_change") = true, py::arg("old_style") = false, py::arg("verbose") = false)
          .def("print_version", &Simulation::print_version);
  }  

}

PYBIND11_MODULE(vm, m)
{
  VMSim::export_Vec(m);
  VMSim::export_VertexProperty(m);
  VMSim::export_HEProperty(m);
  VMSim::export_EdgeProperty(m);
  VMSim::export_FaceProperty(m);
  VMSim::export_Vertex(m);
  VMSim::export_VertexCirculator(m);
  VMSim::export_Edge(m);
  VMSim::export_HalfEdge(m);
  VMSim::export_Face(m);
  VMSim::export_FaceCirculator(m);
  VMSim::export_Mesh(m);
  VMSim::export_System(m);
  VMSim::export_ForceCompute(m);
  VMSim::export_Integrate(m);
  VMSim::export_Topology(m);
  VMSim::export_Dump(m);
  VMSim::export_Simulation(m);
  
  py::register_exception<std::domain_error>(m, "DomainError");
  py::register_exception<std::invalid_argument>(m, "InvalidArgumentError");
  py::register_exception<std::length_error>(m, "LengthError");
}


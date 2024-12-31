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

namespace VMTutorial
{
  void Simulation::run(int steps, bool topological_change, bool old_style, bool verbose)
  {
    if (topological_change) {
      throw runtime_error("Simulation::run - topological_change is not guaranteed to work, this code has changed a lot since that was last used.");
    }
    
    double progress = 0.0;
    if (verbose) { cout << "Simulation::run - Running simulation for " << steps << " steps" << endl; }
    for (int i = sim_step; i < sim_step + steps; i++)
    {
      if (verbose) { cout << "doing step #" << i << endl; }
      if (topological_change)
      {
          _topology.T1();
      }
      
      if (verbose) { cout << "doing integration" << endl; }
      _integ.apply(verbose);
      _sys.set_topology_change(false);

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
          .def("run", &Simulation::run, py::arg("steps"), py::arg("topological_change") = true, py::arg("old_style") = false, py::arg("verbose") = false)
          .def("print_version", &Simulation::print_version);
  }  

}

PYBIND11_MODULE(vm, m)
{
  VMTutorial::export_Vec(m);
  VMTutorial::export_VertexProperty(m);
  VMTutorial::export_HEProperty(m);
  VMTutorial::export_EdgeProperty(m);
  VMTutorial::export_FaceProperty(m);
  VMTutorial::export_Vertex(m);
  VMTutorial::export_VertexCirculator(m);
  VMTutorial::export_Edge(m);
  VMTutorial::export_HalfEdge(m);
  VMTutorial::export_Face(m);
  VMTutorial::export_FaceCirculator(m);
  VMTutorial::export_Mesh(m);
  VMTutorial::export_System(m);
  VMTutorial::export_ForceCompute(m);
  VMTutorial::export_Integrate(m);
  VMTutorial::export_Topology(m);
  VMTutorial::export_Dump(m);
  VMTutorial::export_Simulation(m);
  
  py::register_exception<std::domain_error>(m, "DomainError");
  py::register_exception<std::invalid_argument>(m, "InvalidArgumentError");
  py::register_exception<std::length_error>(m, "LengthError");
}


/*!
 * \file system.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 30-Nov-2023
 * \brief System class 
 */ 

#ifndef __SYSTEM_HPP__
#define __SYSTEM_HPP__

#include <string>
#include <map>
#include <memory>
#include <cmath>

// #include "rng.hpp"
#include "mesh.hpp"
#include "property.hpp"


using std::string;
using std::map;

namespace VMTutorial
{
  typedef Mesh MyMesh;
  

  typedef map<string, double> params_type;         // Used when we set a numerical value to a parameter, e.g., kappa = 1.0
  typedef map<string, Vec> vec_params_type;        // Used when we set a vectorial value to a parameter, e.g., n = Vec(1,0)
  typedef map<string, vector<double>> multi_params_type;   // Used when we need to at least two values to set a parameter, e.g., a parameters is drawn from a random distribution 
  typedef map<string, string> string_params_type;         // Used when we set a string value to a parameter, e.g., force_type = "nematic_vm"
  
  bool operator<(const VertexHandle&, const VertexHandle&);

  class System
  {
    public:

      System(MyMesh& mesh) : _mesh{mesh}, 
                             _mesh_set{false},
                             _topology_changed{true}
                             {
                             }

      // System setup
      void read_input_from_jsonstring(const string& json_contents, bool verbose=false);
      
      void log_debug_stats()
      {
        std::cout << "   System::log_debug_stats:" << std::endl;
        std::cout << "     CURRENT size of _halfedges: " << _mesh.halfedges().size() << std::endl; 
      }
      void set_topology_change(bool flag) { _topology_changed = flag; }

      // System info access 
      MyMesh& mesh()  { return _mesh; }
      const MyMesh& cmesh() const { return _mesh; }
      
      // int& time_step() { return _time_step; }
      // double& simulation_time() { return _simulation_time; }
      bool topology_change() { return _topology_changed;  }
      
    private:
  
      MyMesh &_mesh;
      // int _time_step;
      // double _simulation_time;
      bool _mesh_set;
      bool _topology_changed;  // If true, mesh topology has changed

  };

  void export_VertexProperty(py::module&);
  void export_EdgeProperty(py::module&);
  void export_HEProperty(py::module&);
  void export_FaceProperty(py::module&);
  // void export_Spoke(py::module &);
  void export_Vertex(py::module &);
  void export_VertexCirculator(py::module &);
  void export_Edge(py::module&);
  void export_HalfEdge(py::module&);
  void export_Face(py::module&);
  void export_FaceCirculator(py::module &);
  void export_Mesh(py::module&);
  void export_System(py::module&);

}

#endif

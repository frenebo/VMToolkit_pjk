/*!
 * \file topology.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com 
 * \date 30-Nov-2023
 * \brief Topology class 
*/

#include "topology.hpp"

using std::cout;
using std::endl;

namespace VMTutorial
{
  
  void Topology::T1(bool verbose)
  {
    if (_min_edge_len <= 0 || _new_edge_len <= 0) {
      throw runtime_error("Invalid params - have these been changed from default -1 values? min_edge_len=" + std::to_string(_min_edge_len) + " new_edge_len" + std::to_string(_new_edge_len));
    }
    
    for (auto e : _sys.mesh().edges())
    {
        if (
          // !e.erased &&
          _sys.mesh().len(e) < _min_edge_len
        ) {
          // if (verbose) {
            // cout << " attempting T1 TRANSITION for edge " << e.idx() << "- old len " << _sys.mesh().len(e) << "  ,new len " << _new_edge_len << endl;
          // }
          if (_sys.mesh().T1(e, _new_edge_len, verbose)) {
            cout << "  T1 transition returned true - setting topology change to true" << endl;
            _sys.set_topology_change(true);
          }
        }
    }
  }

  
  void export_Topology(py::module& m)
  {
    py::class_<Topology>(m, "Topology")
      .def(py::init<System&, int>(), py::arg("sys"),  py::arg("seed") = 0)
      .def("set_params", &Topology::set_params)
      ;
  }
  
}

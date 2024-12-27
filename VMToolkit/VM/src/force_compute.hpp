/*!
 * \file force_compute.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 30-Nov-2023
 * \brief ForceCompute class 
 */ 

#ifndef __FORCE_COMPUTE_HPP__
#define __FORCE_COMPUTE_HPP__

#include <exception>
#include <algorithm>
#include <map>
#include <string>
#include <optional>
#include <chrono>

#include "system.hpp"
#include "class_factory.hpp"
#include "force.hpp"
#include "force_area.hpp"
#include "force_perimeter.hpp"
#include "force_const_vertex_propulsion.hpp"
#include "force_efield_on_cell_boundary.hpp"
// #include "force_self_propulsion.hpp"

 
using std::runtime_error;
using std::transform;
using std::map;
using std::string;


namespace VMTutorial
{

  class ForceCompute : public ClassFactory<Force>
  {
    typedef std::chrono::duration<
      std::chrono::high_resolution_clock::rep,
      std::chrono::high_resolution_clock::period
    > timer_duration;
    typedef std::chrono::time_point<std::chrono::high_resolution_clock>  time_pt;
    
    public:

      ForceCompute(System& sys) : _sys{sys}
      { 

      }

      ~ForceCompute() = default; 

      ForceCompute(const ForceCompute&) = delete;
      
      void start_force_compute_timers(bool verbose);
      
      map<string, double> get_force_compute_timers_millis(bool verbose);
      
      std::vector<Vec> compute_all_vertex_forces(bool verbose=false);
      
      double tension(HalfEdge<Property>& he);
      
      void set_global_params(const string& force_id, const params_type& num_params, const map<string,string>& str_params,  bool verbose);

      
      void set_face_params_facewise(const string& force_id, const vector<int>& fids, const vector<params_type>& params, bool verbose);
      
      void set_vertex_params_vertexwise(const string& force_id, const vector<int>& vids, const vector<params_type>& params, bool verbose);

      void add_force(const string& force_id, const string& force_type, bool verbose);

    private:
      System& _sys;
      map<string, timer_duration> _force_timers;
  };

  void export_ForceCompute(py::module&);

}
#endif

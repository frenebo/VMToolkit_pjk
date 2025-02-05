/*!
 * \file force_compute.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 30-Nov-2023
 * \brief ForceCompute class 
 */ 

#ifndef __FORCE_COMPUTE_HPP__
#define __FORCE_COMPUTE_HPP__

#include <map>
#include <string>
#include <chrono>

#include "system.hpp"
#include "class_factory.hpp"
#include "forces/force.hpp"
 

namespace VMSim
{
  using std::map;

  class ForceCompute : public ClassFactory<Force>
  {
    typedef std::chrono::duration<
      std::chrono::high_resolution_clock::rep,
      std::chrono::high_resolution_clock::period
    > timer_duration;
    typedef std::chrono::time_point<std::chrono::high_resolution_clock>  time_pt;
    
    public:

      ForceCompute(const System& sys) : _sys{sys}
      { 

      }

      ~ForceCompute() = default; 

      ForceCompute(const ForceCompute&) = delete;
      
      void start_force_compute_timers(bool verbose);
      
      map<std::string, double> get_force_compute_timers_millis(bool verbose);
      
      void compute_all_vertex_forces(vector<Vec>& res_forces, bool verbose);
      map<string, vector<Vec>> for_py_get_all_vertex_instantaneous_forces(bool verbose);
      
      void set_global_params(
        const std::string& force_id,
        const std::map<string, double>& num_params,
        const std::map<string,string>& str_params,
        const map<string, int>& int_params,
        const map<string, vector<double>> flt_array_params,
        bool verbose
      );
      
      void set_face_params_facewise(const std::string& force_id, const vector<int>& fids, const vector<params_type>& params, bool verbose);
      
      void set_vertex_params_vertexwise(const std::string& force_id, const vector<int>& vids, const vector<params_type>& params, bool verbose);

      void add_force(const std::string& force_id, const std::string& force_type, bool verbose);
      
      void delete_force(const std::string& force_id, bool verbose);

    private:
      const System& _sys;
      std::map<std::string, timer_duration> _force_timers;
  };

  void export_ForceCompute(py::module&);

}
#endif


/*!
 * \file force.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 13-Jun-2017
 * \brief Force class 
 */ 

#ifndef __FORCE_HPP__
#define __FORCE_HPP__

#include <vector>
#include <string>
#include <map>

#include "../system.hpp"

double const SMALL_NUMBER = 1e-6;

using std::vector;
using std::string;
using std::runtime_error;


namespace VMSim
{
  // Force on a vertex
  class Force 
  {
    public:                                                                                                 
      Force(const System& sys) : _sys{sys}
      {
      }
      
      virtual ~Force() { }
      
      // virtual void 
      virtual void compute_all_vertex_forces(vector<Vec>& res, bool verbose)
      {
        throw runtime_error("Unimpllemented Force::compute_al_vertex_forces - has not been overriden in child class");
      }
      
      virtual void set_face_params_facewise(const vector<int>& fids, const vector<params_type>& params, bool verbose)
      {
        throw runtime_error("Unimplemented Force::set_face_params_facewise - has not been overridden in child clas apparently.");
      }
      
      virtual void set_vertex_params_vertexwise(const vector<int>& vids, const vector<params_type>& params, bool verbose)
      {
        throw runtime_error("Unimplemented Force::set_vertex_params_vertexwise - has not been overridden in child class apparently.");
      }
      
      virtual void set_global_params(
        const params_type& num_params,
        const std::map<string,string>& str_params,
        const std::map<string, int>& int_params,
        const map<string, vector<double>> flt_array_params,
        bool verbose
      ) {
        throw runtime_error("Unimplemented Force::set_global_params - apparently has not been overriden in child class.");
      }
      
      virtual void set_global_floatarray_params(const map<string, vector<double>>& params, bool verbose)
      {
        throw runtime_error("Unimplemented Force::set_global_floatarray_params - apparently has not been overriden in child class.");
      }
    
    protected:
      const System& _sys;    // Mes
  };
}

#endif

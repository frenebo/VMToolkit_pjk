
/*!
 * \file force.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 13-Jun-2017
 * \brief Force class 
 */ 

#ifndef __FORCE_HPP__
#define __FORCE_HPP__

#include <vector>
#include <algorithm>
#include <string>
#include <map>
#include <exception>
#include <chrono>

#include "system.hpp"

double const SMALL_NUMBER = 1e-6;

using std::vector;
using std::string;
using std::runtime_error;


namespace VMTutorial
{

  // Force on a vertex
  class Force 
  {
    public:                                                                                                 
      Force(const System& sys) : _sys{sys}
      {
      }
      
      virtual ~Force() { }
        
      // computes force on vertex by a given edge
      virtual Vec compute_he_force(const Vertex<Property>&, const HalfEdge<Property>&, bool verbose=false)
      {
        throw runtime_error("Unimplemented Force::compute_he_force - has not been overridden in child clas apparently.");
      }
      
      virtual void set_face_params_facewise(const vector<int>& fids, const vector<params_type>& params, bool verbose)
      {
        throw runtime_error("Unimplemented Force::set_face_params_facewise - has not been overridden in child clas apparently.");
      }
      
      virtual void set_vertex_params_vertexwise(const vector<int>& vids, const vector<params_type>& params, bool verbose)
      {
        throw runtime_error("Unimplemented Force::set_vertex_params_vertexwise - has not been overridden in child class apparently.");
      }
      
      virtual void set_global_params(const params_type& params,const std::map<string,string>& str_params, bool verbose)
      {
        throw runtime_error("Unimplemented Force::set_global_params - apparently has not been overriden in child class.");
      }
    
    protected:
      const System& _sys;    // Mes
  };
}

#endif

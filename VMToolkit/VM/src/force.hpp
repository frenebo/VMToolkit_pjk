
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
#include "rng.hpp"

double const SMALL_NUMBER = 1e-6;

using std::map;
using std::vector;
using std::string;
using std::runtime_error;
using std::cerr;
using std::exception;
using namespace std::chrono;

namespace VMTutorial
{

  // Force on a vertex
  class Force 
  {
    public:
      // Orphaned model options
      // const param_type& gammaLen, const param_type& l0, const param_type& counter
      // _gammaLen(gammaLen), _l0(l0), _counter(counter)
                                                                                                                                
      Force(const System& sys) : _sys{sys}
      { 

      }
      virtual ~Force() { }
        
      // computes force on vertex by a given edge
      virtual Vec compute_he_force(const Vertex<Property>&, const HalfEdge<Property>&, bool verbose=false)
      {
        throw runtime_error("Unimplemented Force::compute_he_force - has not been overridden in child clas apparently.");
      }
      
      // compute edge tension for handling 4-vertices (and moves that computation out of integrator)
      // takes care of correct type of force law that way
      virtual double tension(const HalfEdge<Property>&)
      {
        throw runtime_error("Unimplemented Force::tension - has not been overridden in child clas apparently.");
      }

      // // computed energy on of a given face
      // virtual double energy(const Face<Property>&) = 0;
      
      // // // set all parameters for a given type
      // virtual void set_params(const string&, const params_type&) = 0;
      
      virtual void set_face_params_facewise(const vector<int>& fids, const vector<params_type>& params, bool verbose)
      {
        throw runtime_error("Unimplemented Force::set_face_params_facewise - has not been overridden in child clas apparently.");
      }
      
      virtual void set_vertex_params_vertexwise(const vector<int>& vids, const vector<params_type>& params, bool verbose)
      {
        throw runtime_error("Unimplemented Force::set_vertex_params_vertexwise - has not been overridden in child class apparently.");
      }
      
      virtual void set_global_params(const params_type& params,const map<string,string>& str_params, bool verbose)
      {
        throw runtime_error("Unimplemented Force::set_global_params - apparently has not been overriden in child class.");
      }
    
    protected:

      const System& _sys;    // Mesh
     
      
  };

  
}

#endif

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

#include "system.hpp"
#include "class_factory.hpp"
#include "force.hpp"
#include "force_area.hpp"
#include "force_perimeter.hpp"
#include "force_const_vertex_propulsion.hpp"
// #include "force_self_propulsion.hpp"

 
using std::runtime_error;
using std::transform;
using std::map;
using std::string;


namespace VMTutorial
{

  class ForceCompute : public ClassFactory<Force>
  {
    public:

      ForceCompute(System& sys) : _sys{sys}
      { 

      }

      ~ForceCompute() = default; 

      ForceCompute(const ForceCompute&) = delete;

      
      // @TODO remove this from stuff - does compute change the force value or not?
      // might accidentally use the wrong one and change it twice at some point...
      Vec compute_v_force(Vertex<Property> &v, bool verbose=false)
      {
        if (verbose)
        {
          cout << "ForceCompute::compute_v_force(vertex) - Computing force on vector... " << v.id << endl;
        }
        Vec tot_force = Vec(0.0,0.0);
        // v.data().force = Vec(0.0,0.0);
        for (auto he : v.circulator()) {
          tot_force += this->compute_he_force(v, he, verbose);
          // v.data().force += this->compute(v, he, verbose);
        }
        if (verbose)
        {
          cout << " TOTAL force on vector " << v.id << ": " << v.data().force.x << ", " << v.data().force.y << endl;
        }
        return tot_force;
      }

      Vec compute_he_force(Vertex<Property> &v, const HalfEdge<Property> &he, bool verbose=false)
      {
        
        if (verbose)
        {
          cout << "  ForceCompute::compute_he_force(vertex, he) - computing force for vertex " << v.id << " half edge " << he.idx() << endl;
        }
        Vec ftot(0,0);
        for (auto& f : this->factory_map) {
          ftot += f.second->compute_he_force(v, he, verbose);
        }
        if (verbose)
        {
          cout << "    TOT force: "<<ftot.x<<","<<ftot.y<<" for vertex " << v.id << " half edge " << he.idx() << endl;
        }
        return ftot;
      }

      double tension(HalfEdge<Property>& he)
      {
        double T = 0.0;
        for (auto& f : this->factory_map)
          T += f.second->tension(he);
        return T;
      }
      
      void set_global_params(const string& force_id, const vector<params_type>&params, bool verbose)
      {
        // throw runtime_error("Unimplemented ForceCompute::set_global_params");
        
        if (this->factory_map.find(force_id) != this->factory_map.end()) {
          this->factory_map[force_id]->set_global_params(params, verbose);
        } else {
          throw runtime_error("ForceCompute::set_global_params: Force id " + force_id + " has not been added yet, could not set its params");
        }
      }

      
      void set_face_params_facewise(const string& force_id, const vector<int>& fids, const vector<params_type>& params, bool verbose)
      {
        if (!(fids.size() == params.size())) {
          throw runtime_error(
            "ForceCompute::set_face_params_facewise: For setting face parameters elementwise, number of fids and fparams should be identical."
          );
        }
        
        if (this->factory_map.find(force_id) != this->factory_map.end()) {
            this->factory_map[force_id]->set_face_params_facewise(fids, params, verbose);
        }
        else {
          throw runtime_error("ForceCompute::set_params_facewise: Force id " + force_id + " has not been added.");
        }
      }
      
      void set_vertex_params_vertexwise(const string& force_id, const vector<int>& vids, const vector<params_type>& params, bool verbose)
      {
        if (!(vids.size() == params.size())) {
          throw runtime_error(
            "ForceCompute::set_vertex_params_vertexwisee: For setting vertex parameters elementwise, number of vids and fparams should by identical."
          );
        }
        
        if (this->factory_map.find(force_id) != this->factory_map.end()) {
          this->factory_map[force_id]->set_vertex_params_vertexwise(vids, params, verbose);
        } else {
          throw runtime_error("ForceCompute::set_vertex_params_vertexwises: Force type " + force_id + " has not been setup in this simulation.");
        }
      }

      void add_force(const string& force_id, const string& force_type, bool verbose)
      {
        if (verbose) {
          cout << "ForceCompute::add_force - Adding force with id '" << force_id << "' and force type '" << force_type << "'" << endl;
        }
        // Check if this force has already been added.
        if (this->factory_map.find(force_id) != this->factory_map.end()) {
          throw runtime_error("ForceCompute::add_force - Cannot add force with id '" + force_id + "' - already has been added to ForceCompute class.");
        }
        
        // if (verbo)
        
        if (force_type == "area") {
          this->add<ForceArea,System&>(force_id, _sys);
        } else if (force_type == "perimeter") {
          this->add<ForcePerimeter,System&>(force_id, _sys);
        } else if (force_type == "const_vertex_propulsion") {
          this->add<ForceConstVertexPropulsion,System&>(force_id, _sys);
        } else  {
          throw runtime_error("Unknown force name : " + force_id + ".");
        }
      }

    private: 

      System& _sys;
 
  };

  void export_ForceCompute(py::module&);

}
#endif

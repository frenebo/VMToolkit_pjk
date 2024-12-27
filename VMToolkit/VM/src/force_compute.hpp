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
      
      void start_force_compute_timers(bool verbose)
      {
        if (verbose)
        {
          cout << "ForceCompute::start_force_compute_timers - clearing and setting timers to zero" << endl;
        }
        _force_timers.clear();
        // _force_timers = map<string, timer_duration>();
        // _force
        
        for ( const auto &it : factory_map ) {
          _force_timers[it.first] = timer_duration::zero();
        }
        if (verbose)
        {
          cout << "  finished initializing timers" << endl;
        }
      }
      
      map<string, double> get_force_compute_timers_millis(bool verbose)
      {
        if (verbose) {
          cout << "ForceCompute::get_force_compute_timers_millis - reading values in timer map" << endl;
        }
        map<string, double> millis_for_forces;
        
        for (const auto &it : _force_timers) {
          string force_id = it.first;
          timer_duration tot_elapsed = it.second;
          
          std::chrono::duration<double, std::milli> ms_double = tot_elapsed;
          
          millis_for_forces[force_id] = ms_double.count();
        }
        
        return millis_for_forces;
      }
      
      // void compute_and_ap
      void compute_and_set_all_vertex_forces(bool verbose=false)
      {
        if (verbose) {
          cout << "ForceCompute::compute_and_apply_vertex_force - computing all forces" << endl;
        }
        
          
        std::vector<Vec> v_forces(_sys.mesh().vertices().size(), Vec(0.0,0.0));
        
        bool do_time_computations = (_force_timers.size() > 1);
        
          // Vector forces start at zero, we add force as w iterate through forces, vectors, halfedges
        for (auto& fit : this->factory_map)
        {
          std::string force_id = fit.first;
          if (verbose)
          {
            cout << "  computing force '" << force_id << "'" << endl;
          }
          
          auto& force = fit.second;
          
          // 
          time_pt t1;
          time_pt t2;
          bool do_time_this_force = (do_time_computations && _force_timers.count(force_id));
          
          if (do_time_this_force) {
            t1 = std::chrono::high_resolution_clock::now();
          }
          
          for (auto& vertex : _sys.mesh().vertices())
          {
            if (vertex.erased) {
              continue;
            }
            
            for (auto& he : vertex.circulator()) {
              if (verbose)
              {
                cout << "    computing force on vertex " << vertex.id << " by halfedge " << he.idx() << endl;
              }
              
              v_forces.at(vertex.id) += force->compute_he_force(vertex, he, verbose);
            }
          }
          
          if (do_time_this_force) {
            t2 = std::chrono::high_resolution_clock::now();
            auto elapsed = t2 - t1;
            
            _force_timers.at(force_id) += elapsed;
          }
        }
        
        size_t vid = 0;
        for (const Vec& vforce : v_forces) {
          _sys.mesh().vertices().at(vid).data().force = vforce;
          
          vid++;
        }
      }
      
      // // @TODO remove this from stuff - does compute change the force value or not?
      // // might accidentally use the wrong one and change it twice at some point...
      // Vec compute_v_force(Vertex<Property> &v, bool verbose=false)
      // {
      //   if (verbose)
      //   {
      //     cout << "ForceCompute::compute_v_force(vertex) - Computing force on vector... " << v.id << endl;
      //   }
      //   Vec tot_force = Vec(0.0,0.0);
      //   // v.data().force = Vec(0.0,0.0);
      //   for (auto he : v.circulator()) {
      //     tot_force += this->compute_he_force(v, he, verbose);
      //     // v.data().force += this->compute(v, he, verbose);
      //   }
      //   if (verbose)
      //   {
      //     cout << " TOTAL force on vector " << v.id << ": " << v.data().force.x << ", " << v.data().force.y << endl;
      //   }
      //   return tot_force;
      // }

      // Vec compute_he_force(Vertex<Property> &v, const HalfEdge<Property> &he, bool verbose=false)
      // {
        
      //   if (verbose)
      //   {
      //     cout << "  ForceCompute::compute_he_force(vertex, he) - computing force for vertex " << v.id << " half edge " << he.idx() << endl;
      //   }
      //   Vec ftot(0,0);
      //   for (auto& f : this->factory_map) {
      //     ftot += f.second->compute_he_force(v, he, verbose);
      //   }
      //   if (verbose)
      //   {
      //     cout << "    TOT force: "<<ftot.x<<","<<ftot.y<<" for vertex " << v.id << " half edge " << he.idx() << endl;
      //   }
        
      //   return ftot;
      // }

      double tension(HalfEdge<Property>& he)
      {
        double T = 0.0;
        for (auto& f : this->factory_map)
          T += f.second->tension(he);
        return T;
      }
      
      void set_global_params(const string& force_id, const params_type& num_params, const map<string,string>& str_params,  bool verbose)
      {
        // throw runtime_error("Unimplemented ForceCompute::set_global_params");
        
        if (this->factory_map.find(force_id) != this->factory_map.end()) {
          this->factory_map[force_id]->set_global_params(num_params,str_params, verbose);
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
        
        if (force_type == "area") {
          this->add<ForceArea,System&>(force_id, _sys);
        } else if (force_type == "perimeter") {
          this->add<ForcePerimeter,System&>(force_id, _sys);
        } else if (force_type == "const_vertex_propulsion") {
          this->add<ForceConstVertexPropulsion,System&>(force_id, _sys);
        } else if (force_type == "force_efield_on_cell_boundary") {
          this->add<ForceEFieldOnCellBoundary, System&>(force_id, _sys);
        } else  {
          throw runtime_error("Unknown force name : " + force_id + ".");
        }
      }

    private:
      System& _sys;
      map<string, timer_duration> _force_timers;
  };

  void export_ForceCompute(py::module&);

}
#endif

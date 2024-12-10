/*!
 * \file integrator_brownian.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 30-Nov-2023
 * \brief IntegratorBrownian class 
*/

#ifndef __INTEGRATOR_BROWNIAN_HPP__
#define __INTEGRATOR_BROWNIAN_HPP__

#include "constrainer.hpp"
#include "constraint_none.hpp"
#include "constraint_fixed.hpp"


#include "integrator.hpp"

#include <chrono>
#include <utility>
#include <map>
#include <string> 
#include <memory>

using namespace std::chrono;
using std::map;
using std::make_unique;

namespace VMTutorial
{

  using ConstrainerType = unique_ptr<Constrainer>;

  class IntegratorBrownian : public Integrator 
  {

    public:

      IntegratorBrownian(System& sys, ForceCompute& fc, int seed) : Integrator{sys, fc, seed},
                                                                    _T{0.0},
                                                                    _gamma{1.0},
                                                                    _Dr{0.0},
                                                                    _update_n{false}
                                                                    
      { 
        map<string,int>& vert_types = _sys.vert_types();
        map<string,int>& cell_types = _sys.cell_types();
        for (int i = 0; i < vert_types.size(); i++) {
          _constant_force.push_back(Vec(0.0,0.0));
        }

        _constrainer = make_unique<Constrainer>();
        _constrainer->add<ConstraintNone>("none");
        _constrainer->add<ConstraintFixed>("fixed");
      }

      void step(bool verbose) override;
      void set_params(const params_type& params) override 
      { 
        for (auto& p : params)
        {
          if (p.first == "T")
            _T = p.second;
          if (p.first == "gamma")
            _gamma = p.second;
          if (p.first == "Dr")
            _Dr = p.second;
        }
      };
      void set_type_params(const string& type, const params_type& params) override { }
      void set_string_params(const string_params_type& params) override { }
      void set_external_force(const string& vtype, const Vec& f) override 
      { 
        _constant_force[_sys.vert_types()[vtype]] = f;
      }
      void set_external_forces_by_vertex(const vector<int>& vids, const vector<Vec>& forces) override
      {
        if (vids.size() != forces.size()) {
          throw runtime_error("Number of vertex indices and forces do not match");
        }
        throw runtime_error("This doesn't work the way I thought it did... this sets forces by type");
        // if ()
        
        // for (size_t i = 0; i < vids.size(); ++i)
        // {
        //   if (i >= vids.size()) {
        //     throw runtime_error("there are not enough vids in here. Count is " + std::to_string(vids.size()));
        //   }
        //   if (i >= forces.size()) {
        //     throw runtime_error("there are not enough forces in here. Count is " + std::to_string(forces.size()));
        //   }
        //   std::cout << "Setting force of vertex " << vids[i] << " to " << forces[i].x << ", " << forces[i].y << std::endl;
          
        //   _constant_force.at(vids[i]) = forces[i];
        // }
        // std::cout 
        // throw runtime_error("Not implemented");
      }
      
      void set_flag(const string& flag) override 
      {  
        if (flag == "update_n")
          _update_n = true;
        else
          throw runtime_error("Brownian integrator: Unknown flag : " + flag + ".");
      }

      
    private:

      ConstrainerType _constrainer; // Apply various constraints
      double _T;                 // temperature
      double _gamma;             // friction 
      vector<Vec> _constant_force;
      double _Dr;                // rotational diffusion constant
      bool _update_n;            // If true, update direction of the cell self-propulsion direction
      
    };

}

#endif


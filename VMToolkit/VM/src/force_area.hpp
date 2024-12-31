/*!
 * \file force_area.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 30-Nov-2023
 * \brief ForceArea class 
 */ 

#ifndef __FORCE_AREA_HPP__
#define __FORCE_AREA_HPP__

#include "force.hpp"

#include <vector>
#include <string>


namespace VMTutorial
{
  // Force on a vertex
  class ForceArea : public Force
  {
    public:
                                                                                                                                
      ForceArea(const System& sys) : Force(sys) 
      { 
      }
      virtual ~ForceArea() { }
        
      void compute_all_vertex_forces(vector<Vec>& res, bool verbose) override;
      
      void set_face_params_facewise(const std::vector<int>& fids, const std::vector<params_type>& params, bool verbose) override;
    private:
      Vec compute_he_force(const Vertex&, const HalfEdge&, bool verbose);
      bool _enabled_for_faceidx(int fid, bool verbose);
      
      void _clear_compute_cache(bool verbose);
      void _cache_mesh_computations(bool verbose);
      
      std::vector<bool> _force_enabled_mask_by_cell_index;
      std::vector<double> _kappa_params;
      std::vector<double> _A0_params;
      
      std::vector<double> _cached_face_areas;
      std::vector<bool> _cached_face_enabled;
  };

  
}

#endif

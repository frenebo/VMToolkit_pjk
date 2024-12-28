/*!
 * \file force_area.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 30-Nov-2023
 * \brief ForceArea class 
 */ 

#ifndef __FORCE_AREA_HPP__
#define __FORCE_AREA_HPP__

#include "force.hpp"


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
        
      // computes force on vertex by a given edge
      Vec compute_he_force(const Vertex<Property>&, const HalfEdge<Property>&, bool verbose) override;  
      
      
      void set_face_params_facewise(const vector<int>& fids, const vector<params_type>& params, bool verbose) override;

      
      bool enabled_for_faceidx(int fid, bool verbose);
    
    private:
      vector<bool> _force_enabled_mask_by_cell_index;
      vector<double> _kappa;
      vector<double> _A0;
  };

  
}

#endif

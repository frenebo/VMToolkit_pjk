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
                                                                                                                                
      ForceArea(System& sys) : Force(sys) 
      { 
      }
      virtual ~ForceArea() { }
        
      // computes force on vertex by a given edge
      Vec compute(const Vertex<Property>&, const HalfEdge<Property>&) override;  
      
      // Are term does not generate tension
      double tension(const HalfEdge<Property>& he) override
      {
        return 0.0;
      }


      // Energy calculation 
      double energy(const Face<Property>&) override;
      
      
      void set_face_params_facewise(const vector<int>& fids, const vector<params_type>& params) override {
        // Resize the internal parameter arrays if needed...
        int max_fid = 0;
        for (const auto& fid : fids) {
          if (fid > max_fid) max_fid = fid;
        }
        
        if (max_fid >= _kappa.size()) {
          _kappa.resize(max_fid + 1, 0.0);
        }
        
        for (size_t i = 0; i < fids.size(); ++i) {
          int fid = fids[i];
          const auto& fparam = params[i];
 
          _kappa.at(fid) = fparam.at("kappa");
        }
      }


      void set_flag(const string& flag) override   {
        throw runtime_error("No flags implemented for ForceArea");
      }
    
    private:

      vector<double> _kappa;
      
      
  };

  
}

#endif

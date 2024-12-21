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
      Vec compute_he_force(const Vertex<Property>&, const HalfEdge<Property>&, bool verbose=false) override;  
      
      // Are term does not generate tension
      double tension(const HalfEdge<Property>& he) override
      {
        return 0.0;
      }


      // // Energy calculation 
      // double energy(const Face<Property>&) override;
      
      
      void set_face_params_facewise(const vector<int>& fids, const vector<params_type>& params, bool verbose) override {
        // Resize the internal parameter arrays if needed...
        int max_fid = 0;
        for (const auto& fid : fids) {
          if (fid > max_fid) max_fid = fid;
        }
        
        if (max_fid >= _force_enabled_mask_by_cell_index.size()) {
          _force_enabled_mask_by_cell_index.resize(max_fid + 1, false);
          _kappa.resize(max_fid + 1, 0.0);
          _A0.resize(max_fid + 1, 0.0);
        }
        
        for (size_t i = 0; i < fids.size(); ++i) {
          int fid = fids[i];
          const auto& fparam = params[i];
          
          _force_enabled_mask_by_cell_index.at(fid) = true;
          _kappa.at(fid) = fparam.at("kappa");
          _A0.at(fid) = fparam.at("A0");
          
          if (verbose) {
            cout << "Setting area force for face '" << fid << "':  kappa=" << fparam.at("kappa") << " A0=" << fparam.at("A0") << endl;
          }
        }
      }


      // void set_flag(const string& flag) override   {
      //   throw runtime_error("No flags implemented for ForceArea");
      // }
      
      bool enabled_for_faceidx(int fid, bool verbose=false) {
        // return 
        if (fid >= _force_enabled_mask_by_cell_index.size()) return false;
        
        bool is_enabled =  _force_enabled_mask_by_cell_index[fid];
        
        if (verbose)
        {
          if (is_enabled)
          {
            cout << "          - Confirmed that area force is enabled for fid=" << fid << endl;
          }
          else
          {
            cout << "          - Area force NOT enabled for fid=" << fid << endl;
          }
        }
      
        return is_enabled;
      }
    
    private:
      vector<bool> _force_enabled_mask_by_cell_index;
      vector<double> _kappa;
      vector<double> _A0;
      
      
  };

  
}

#endif

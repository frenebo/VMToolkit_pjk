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
      
      // // set all parameters for a given type
      // void set_params(const string& cell_type, const params_type& params) override
      // {
      //   for (auto p : params)
      //     if (p.first != "kappa")
      //       throw runtime_error("Unknown parameter "+p.first+".");
            
      //   if (params.find("kappa") == params.end())
      //     throw runtime_error("Area force requires parameter kappa.");

      //   try 
      //   {
      //     if (_sys.cell_types().find(cell_type) == _sys.cell_types().end()) {
      //       cout << "Cell types are :";
      //       for(auto iter = _sys.cell_types().begin(); iter != _sys.cell_types().end(); ++iter)
      //       {
      //         // asdfasdf;
      //         auto k =  iter->first;
      //         auto v = iter->second;
      //         cout << v << " ";
      //       //ignore value
      //       //Value v = iter->second;
      //       }
      //       cout << endl;
      //       throw runtime_error("Force area: Cell type " + cell_type + " is not defined.");
      //     }
      //     if (_kappa.size() < _sys.get_num_cell_types())
      //       _kappa.resize(_sys.get_num_cell_types(), 0.0);
      //     int ct = _sys.cell_types()[cell_type];
      //     _kappa[ct] = params.at("kappa");
      //   } 
      //   catch(const exception& e)
      //   {
      //     cerr << "Problem with setting area force parameters. Exception: " << e.what() << '\n';
      //     throw;
      //   }
      // }
      
      // void _set_params_for_face(int fid, params_type& params) {
      //   // _kappa[]
      //   // _kappa.at(fid)
      // }
      
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
          const params_type& fparam = params[i];
 
          // _set_params_for_face(fid, fparam);
          _kappa.at(fid) = fparam.at("kappa");
        }
        // throw runtime_error("unimplemented");
      }


      void set_flag(const string& flag) override   {
        throw runtime_error("No flags implemented for ForceArea");
      }
    
    private:

      vector<double> _kappa; 
      
      
  };

  
}

#endif

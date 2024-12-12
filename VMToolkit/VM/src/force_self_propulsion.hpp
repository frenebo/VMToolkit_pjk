// /*!
//  * \file force_self_propulsion.hpp
//  * \author Rastko Sknepnek, sknepnek@gmail.com
//  * \date 30-Nov-2023
//  * \brief ForceSelfPropulsion class 
//  */ 

// #ifndef __FORCE_SELF_PROPULSION_HPP__
// #define __FORCE_SELF_PROPULSION_HPP__

// #include "force.hpp"


// namespace VMTutorial
// {

//   // Force on a vertex
//   class ForceSelfPropulsion : public Force
//   {
//     public:
                                                                                                                                
//       ForceSelfPropulsion(System& sys) : Force{sys} 
//       { 
//         // @TODO this should automatically resize when an index bigger that seen before is added
//         _v0.resize(_sys.cell_types().size(), 0.0);
//         _n.resize(_sys.cell_types().size(), Vec(0,0));
//       }

//       virtual ~ForceSelfPropulsion() { }
        
//       // computes force on vertex by a given edge
//       Vec compute(const Vertex<Property>&, const HalfEdge<Property>&) override;  
      
      
//       double tension(const HalfEdge<Property>& he) override
//       {
//         return 0.0;
//       }

      
//       // Energy of self propulsion is not defined, so we return 0
//       double energy(const Face<Property>& f) override
//       {
//         return 0.0;   
//       }
      
//       // set all parameters for a given type
//       void set_params(const string& cell_type, const params_type& params) override
//       {
//         for (auto p : params)
//           if (p.first != "v0")
//             throw runtime_error("Unknown parameter " + p.first + ".");
            
//         if (params.find("v0") == params.end())
//           throw runtime_error("SelfPropulsion force requires parameter v0.");

//         try 
//         {
//           if (_sys.cell_types().find(cell_type) == _sys.cell_types().end())
//             throw runtime_error("Force self-propulsion: Cell type " + cell_type + " is not defined.");
//           if (_v0.size() < _sys.get_num_cell_types())
//             _v0.resize(_sys.get_num_cell_types(), 0.0);
//           int ct = _sys.cell_types()[cell_type];
//           _v0[ct] = params.at("v0");
//         } 
//         catch(const exception& e)
//         {
//           cerr << "Problem with setting self-propulsion force parameters. Exception: " << e.what() << '\n';
//           throw;
//         }
//       }
      
//       void set_face_params_facewise(const vector<int>& fids, const vector<params_type>& params) override {
//         // Expand the params if necessary
//         int max_fid = 0;
//         for (const auto& fid : fids) {
//         if (fid > max_fid) max_fid = fid;
//         }
        
//         if (max_fid >= _v0.size()) {
//           _v0.resize(max_fid + 1, 0.0);
//         }
//         if (max_fid >= _n.size()) {
//           _n.resize(max_fid + 1, new Vec(0.0,0.0));)
//         }
        
//         // Set the values
//         for (size_t i =0; i < fids.size(); i++) {
//           int fid = fids[i];
//           params_type& fparam = params[i];
          
//           _v0.at(fid) = fparam.at("v0");
//         }
//       }

//       // set all vector-valued parameters for a given type
//       void OLD_set_vec_params(const string& cell_type, const vec_params_type& params) override 
//       { 
//         for (auto p : params)
//           if (p.first != "n")
//             throw runtime_error("Unknown parameter " + p.first + ".");

//         if (params.find("n") == params.end())
//           throw runtime_error("SelfPropulsion force requires vector parameter n.");

//         try 
//         {
//           int ct = _sys.cell_types()[cell_type];
//           if (ct + 1 > _n.size())
//             _n.push_back(params.at("n").unit());
//           else
//             _n[ct] = params.at("n").unit();
//         } 
//         catch(const exception& e)
//         {
//           cerr << "Problem with setting self-propulsion force parameters. Exception: " << e.what() << '\n';
//           throw e;
//         }

//       };

//       void set_flag(const string& flag) override   {
//         throw runtime_error("No flags implemented for ForceSelfPropulsion");
//       }
    
//     private:

//       vector<double> _v0; 
//       // vector<Vec> _n; 
      
      
//   };

  
// }

// #endif

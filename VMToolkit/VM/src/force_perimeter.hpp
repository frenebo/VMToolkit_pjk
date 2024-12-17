/*!
 * \file force_perimeter.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 30-Nov-2023
 * \brief ForcePerimeter class
 */

#ifndef __FORCE_PERIMETER_HPP__
#define __FORCE_PERIMETER_HPP__

#include "force.hpp"

namespace VMTutorial
{

	// Force on a vertex
	class ForcePerimeter : public Force
	{
	public:
		ForcePerimeter(System &sys) : Force{sys}
		{
			// @TODO this should automatically resize when an index bigger that seen before is added
			
			// _gamma.resize(_sys.cell_types().size(), 0.0);
			// _lambda.resize(_sys.cell_types().size(), 0.0);
		}
		virtual ~ForcePerimeter() {}

		// computes force on vertex by a given edge
		Vec compute(const Vertex<Property> &, const HalfEdge<Property> &, bool verbose=false) override;

		double tension(const HalfEdge<Property> &) override;

		// // Energy calculation
		// double energy(const Face<Property> &) override;

		// set all parameters for a given type
		// void set_params(const string &cell_type, const params_type &params) override
		// {
		// 	for (auto p : params)
		// 		if (p.first != "gamma" && p.first != "lambda")
		// 			throw runtime_error("Unknown parameter " + p.first + ".");

		// 	if (params.find("gamma") == params.end())
		// 		throw runtime_error("Perimeter force requires parameter gamma.");
		// 	if (params.find("lambda") == params.end())
		// 		throw runtime_error("Perimeter force requires parameter lambda.");

		// 	try
		// 	{
		// 		if (_sys.cell_types().find(cell_type) == _sys.cell_types().end())
		// 			throw runtime_error("Force perimeter: Cell type " + cell_type + " is not defined.");
		// 		if (_gamma.size() < _sys.get_num_cell_types())
		// 			_gamma.resize(_sys.get_num_cell_types(), 0.0);
		// 		if (_lambda.size() < _sys.get_num_cell_types())
		// 			_lambda.resize(_sys.get_num_cell_types(), 0.0);
		// 		int ct = _sys.cell_types()[cell_type];
		// 		_gamma.at(ct) = params.at("gamma");
		// 		_lambda.at(ct) = params.at("lambda");
		// 	}
		// 	catch (const exception &e)
		// 	{
		// 		cerr << "Problem with setting perimeter force parameters. Exception: " << e.what() << '\n';
		// 		throw;
		// 	}
		// }

		void set_face_params_facewise(const vector<int>& fids, const vector<params_type>& params, bool verbose) override {
			// Expand the params if necessary
			int max_fid = 0;
			for (const auto& fid : fids) {
			if (fid > max_fid) max_fid = fid;
			}
			
			if (max_fid >= _force_enabled_mask_by_cell_index.size()) {
				_force_enabled_mask_by_cell_index.resize(max_fid + 1, false);
				_gamma.resize(max_fid + 1, 0.0);
				_lambda.resize(max_fid + 1, 0.0);
			}
			
			// Set the param values in the vectors
			for (size_t i = 0; i < fids.size(); ++i) {
				int fid = fids[i];
				const params_type& fparam = params.at(i);
				
				_force_enabled_mask_by_cell_index.at(fid) = true;
				_gamma.at(fid) = fparam.at("gamma");
				_lambda.at(fid) = fparam.at("lambda");
			}
		}
      
		bool enabled_for_faceidx(int fid, bool verbose) {
			if (fid >= _force_enabled_mask_by_cell_index.size()) return false;
			
			bool is_enabled = _force_enabled_mask_by_cell_index[fid];
			
			if (verbose)
			{
				if (is_enabled) {
					cout << "Confirmed that perimeter force is on for fid=" << fid << endl;
				}
			}
			
			return is_enabled;
		}
    
    private:
		vector<bool> _force_enabled_mask_by_cell_index;
		vector<double> _gamma;
		vector<double> _lambda;
	};

}

#endif

/*!
 * \file force_area.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 19-May-2019
 * \brief ForceArea class
 */

#include "force_area.hpp"

using std::cout;
using std::endl;

namespace VMTutorial
{
	Vec ForceArea::compute_he_force(const Vertex<Property> &v, const HalfEdge<Property> &he, bool verbose)
	{
		if (verbose)
		{
			cout << "        ForceArea::compute - computing force for vertex " << v.id << ", halfedge idx " << he.idx()  << endl;
		}
		
		
		Vec l = he.to()->data().r - v.data().r;					// vector along the junction pointing away from the vertex
		const Face<Property> &f = *(he.face());			// cell to the right of the half edge
		const Face<Property> &fp = *(he.pair()->face()); // pair cell (opposite side of the same junction)
		
		bool enabled_for_f = enabled_for_faceidx(f.id, verbose);
		bool enabled_for_fp = enabled_for_faceidx(fp.id, verbose);
		
		
		// If neither adjoining face have this force enabled, then this edge will exert no areaforce on any vertex
		if ((!enabled_for_f) && (!enabled_for_fp)) {
			return Vec(0.0,0.0);
		}
		
		double A1 = _sys.cmesh().area(f);
		double A2 = _sys.cmesh().area(fp);
		
		double kappa_1 = 0.0;
		double A0_1 = 0.0;
		if (enabled_for_f) {
			kappa_1 = (f.outer) ? 0.0 : _kappa.at(f.id);
			A0_1 = _A0.at(f.id);
		}
		
		double kappa_2 = 0.0;
		double A0_2 = 0.0;
		if (enabled_for_fp) {
			kappa_2 = (fp.outer) ? 0.0 : _kappa.at(fp.id);
			A0_2 = _A0.at(fp.id);
		}

		Vec farea_vec = 0.5 * (kappa_1 * (A1 - A0_1) - kappa_2 * (A2 - A0_2)) * l.ez_cross_v();
		
		if (verbose)
		{
			cout << "            A1=" << A1 <<", A2=" << A2 << ", A0_1=" << A0_1 << ", A0_2=" << A0_2 << endl;
			cout << "            area force is (" << farea_vec.x << ", " << farea_vec.y << ")" << endl;
		}

		return farea_vec;
	}

	void ForceArea::set_face_params_facewise(const vector<int>& fids, const vector<params_type>& params, bool verbose)
	{
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
	
	bool ForceArea::enabled_for_faceidx(int fid, bool verbose)
	{
		if (fid >= _force_enabled_mask_by_cell_index.size())
		{
			return false;
		}
		
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
}

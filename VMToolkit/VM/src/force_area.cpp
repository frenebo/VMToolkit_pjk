/*!
 * \file force_area.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 19-May-2019
 * \brief ForceArea class
 */

#include "force_area.hpp"



using std::vector;
using std::cout;
using std::endl;

namespace VMTutorial
{
	void ForceArea::compute_all_vertex_forces(vector<Vec>& res, bool verbose)
	{
		if (verbose) {
			cout << "   ForceArea::compute_all_vertex_forces - starting" << endl;
		}
		_clear_compute_cache(verbose);
		_cache_mesh_computations(verbose);
		
		size_t n_vertices = _sys.cmesh().cvertices().size();	
		res.resize(n_vertices, Vec(0.0,0.0));
			
		for (auto& vertex : _sys.cmesh().cvertices())
		{
			Vec v_force(0.0,0.0);
			for (auto& he : vertex.circulator()) {
				if (verbose) {
					cout << "    computing force on vertex " << vertex.id << " by halfedge " << he.idx() << endl;
				}
				
				v_force += compute_he_force(vertex, he, verbose);
			}
			res.at(vertex.id) = v_force;
		}
	}
	
	void ForceArea::_clear_compute_cache(bool verbose)
	{
		_cached_face_areas.clear();
		_cached_face_enabled.clear();
	}
	
	void ForceArea::_cache_mesh_computations(bool verbose)
	{
		size_t n_faces = _sys.cmesh().cfaces().size();	
		
		_cached_face_areas.resize(n_faces, 0.0);
		_cached_face_enabled.resize(n_faces, false);
		
		size_t fid = 0;
		for (const auto& face : _sys.cmesh().cfaces()) {
			double A = _sys.cmesh().area(face);
			_cached_face_areas.at(fid) = A;
			
			bool f_enabled = _enabled_for_faceidx(fid, verbose);
			_cached_face_enabled.at(fid) = f_enabled;
			
			fid++;
		}
	}
	
	Vec ForceArea::compute_he_force(const Vertex &v, const HalfEdge &he, bool verbose)
	{
		if (verbose)
		{
			cout << "        ForceArea::compute - computing force for vertex " << v.id << ", halfedge idx " << he.idx()  << endl;
		}
		
		
		Vec l = he.to()->data().r - v.data().r;					// vector along the junction pointing away from the vertex
		
		const Face &f = *(he.face());			// cell to the right of the half edge
		const Face &fp = *(he.pair()->face()); // pair cell (opposite side of the same junction)
		
		bool enabled_for_f = _cached_face_enabled[f.id];
		bool enabled_for_fp = _cached_face_enabled[fp.id];
		
		
		// If neither adjoining face have this force enabled, then this edge will exert no areaforce on any vertex
		if ((!enabled_for_f) && (!enabled_for_fp)) {
			return Vec(0.0,0.0);
		}
		
		double A1 = _cached_face_areas[f.id];
		double A2 = _cached_face_areas[fp.id];
		
		double kappa_1 = 0.0;
		double A0_1 = 0.0;
		if (enabled_for_f) {
			kappa_1 = (f.outer) ? 0.0 : _kappa_params.at(f.id);
			A0_1 = _A0_params.at(f.id);
		}
		
		double kappa_2 = 0.0;
		double A0_2 = 0.0;
		if (enabled_for_fp) {
			kappa_2 = (fp.outer) ? 0.0 : _kappa_params.at(fp.id);
			A0_2 = _A0_params.at(fp.id);
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
			_kappa_params.resize(max_fid + 1, 0.0);
			_A0_params.resize(max_fid + 1, 0.0);
		}
		
		for (size_t i = 0; i < fids.size(); ++i) {
			int fid = fids[i];
			const auto& fparam = params[i];
			
			_force_enabled_mask_by_cell_index.at(fid) = true;
			_kappa_params.at(fid) = fparam.at("kappa");
			_A0_params.at(fid) = fparam.at("A0");
			
			if (verbose) {
				cout << "Setting area force for face '" << fid << "':  kappa=" << fparam.at("kappa") << " A0=" << fparam.at("A0") << endl;
			}
		}
	}
	
	bool ForceArea::_enabled_for_faceidx(int fid, bool verbose)
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

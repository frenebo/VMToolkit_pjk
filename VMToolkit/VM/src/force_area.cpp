/*!
 * \file force_area.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 19-May-2019
 * \brief ForceArea class
 */

#include "force_area.hpp"

namespace VMTutorial
{
	Vec ForceArea::compute(const Vertex<Property> &v, const HalfEdge<Property> &he, bool verbose)
	{
		if (verbose)
		{
			cout << "ForceArea::compute - computing force for vertex " << v.id << ", halfedge idx " << he.idx()  << endl;
		}
		
		
		// if (_force_enabled_mask_by_cell_index.)
		Vec l = he.to()->r - v.r;					// vector along the junction pointing away from the vertex
		const Face<Property> &f = *(he.face());			// cell to the right of the half edge
		const Face<Property> &fp = *(he.pair()->face()); // pair cell (opposite side of the same junction)
		double A1 = _sys.mesh().area(f);
		double A2 = _sys.mesh().area(fp);
		// double A0_1 = f.data().A0;
		
		double kappa_1 = 0.0;
		double A0_1 = 0.0;
		if (enabled_for_faceidx(f.id, verbose)) {
			kappa_1 = (f.outer) ? 0.0 : _kappa.at(f.id);
			A0_1 = _A0.at(f.id);
		}
		
		double kappa_2 = 0.0;
		double A0_2 = 0.0;
		if (enabled_for_faceidx(fp.id, verbose)) {
			kappa_2 = (fp.outer) ? 0.0 : _kappa.at(fp.id);
			A0_2 = _A0.at(fp.id);
		}

		Vec farea_vec = 0.5 * (kappa_1 * (A1 - A0_1) - kappa_2 * (A2 - A0_2)) * l.ez_cross_v();
		
		if (verbose)
		{
			cout << "    A1=" << A1 <<", A2=" << A2 << ", A0_1=" << A0_1 << ", A0_2=" << A0_2 << endl;
			cout << "    area force is (" << farea_vec.x << ", " << farea_vec.y << ")" << endl;
		}

		return farea_vec;
	}

	// double ForceArea::energy(const Face<Property> &f)
	// {
	// 	double A = _sys.mesh().area(f);
		
	// 	if (f.outer || f.erased)
	// 		return 0.0;

	// 	double A0;
	// 	double kappa = 0.0;
	// 	if (enabled_for_faceidx(f.id)) {
	// 		A0  = _A0.at(f.id);
	// 		kappa = _kappa.at(f.id);
	// 	}

	// 	double dA = A - A0;
	// 	return 0.5 * kappa * dA * dA;
	// }
}

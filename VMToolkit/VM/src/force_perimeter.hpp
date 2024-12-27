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
		}
		virtual ~ForcePerimeter() {}

		// computes force on vertex by a given edge
		Vec compute_he_force(const Vertex<Property> &, const HalfEdge<Property> &, bool verbose=false) override;

		void set_face_params_facewise(const vector<int>& fids, const vector<params_type>& params, bool verbose) override;
      
		bool enabled_for_faceidx(int fid, bool verbose);
    
    private:
		vector<bool> _force_enabled_mask_by_cell_index;
		vector<double> _gamma;
		vector<double> _lambda;
	};

}

#endif

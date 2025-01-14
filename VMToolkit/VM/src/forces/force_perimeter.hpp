/*!
 * \file force_perimeter.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 30-Nov-2023
 * \brief ForcePerimeter class
 */

#ifndef __FORCE_PERIMETER_HPP__
#define __FORCE_PERIMETER_HPP__

#include "force.hpp"

namespace VMSim
{

	// Force on a vertex
	class ForcePerimeter : public Force
	{
	public:
		ForcePerimeter(const System &sys) : Force{sys}
		{
		}
		virtual ~ForcePerimeter() {}

		// computes force on vertex by a given edge
		void compute_all_vertex_forces(vector<Vec>& res, bool verbose) override;
		void set_face_params_facewise(const vector<int>& fids, const vector<params_type>& params, bool verbose) override;
      
		bool enabled_for_faceidx(int fid, bool verbose);
    
    private:
		Vec _compute_he_force(const Vertex &, const HalfEdge &, bool verbose=false);
		
		vector<bool> _force_enabled_mask_by_cell_index;
		vector<double> _gamma;
		vector<double> _lambda;
	};

}

#endif

/*!
 * \file integrator_adaptive_runge_kutta.cpp
 * \author Paul Kreymborg, frenebo@gmail.com
 * \date 31-Dec-2024
 * \brief IntegratorAdaptiveRungeKutta class 
*/

#include "integrator_adaptive_runge_kutta.hpp"

using std::cout;
using std::endl;

namespace VMTutorial
{
  vector<Vec> IntegratorAdaptiveRungeKutta::_calc_instantaneous_velocities(vector<Vec>& vertex_vels_out, bool verbose) const
  {
    double mu = 1.0 / _gamma_param;    // mobility 
    if (verbose) {
      cout << "  calculating instantaneous velocities with mu=1.0/_gamma_param=" << mu << endl;
    }
    
    _force_compute.compute_all_vertex_forces(vertex_vels_out, verbose);
    
    for (auto& vel : vertex_vels_out)
    {
      vel.x *= mu;
      vel.y *= mu;
    }
    
    return vertex_vels_out;
  }
  
  const vector<double>& IntegratorAdaptiveRungeKutta::_getkvec(int kidx) const {
    if (kidx == 0) return _k0;
    if (kidx == 1) return _k1;
    if (kidx == 2) return _k2;
    if (kidx == 3) return _k3;
    if (kidx == 4) return _k4;
    if (kidx == 5) return _k5;
    if (kidx == 6) return _k6;
    
    throw runtime_error("invalid k idx" + std::to_string(kidx));
  }
  
  void IntegratorAdaptiveRungeKutta::calculate_kvec(int kidx, double hstep, const vector<Vec> & init_positions, bool verbose) {
    // currently unused
    // double dt = hstep * _butcher_tableau_timeoffsets.at(kidx);
    
    const vector<double>& k_compute_coeffs = _butcher_tableau_kcalc_coeffs.at(kidx);
    
    int n_vertices = _sys.mesh().vertices().size();
    
    
    
    // vector<Vec> starting_positions
    for (int vidx = 0; vidx < n_vertices; vidx++) {
      Vec starting_pos = init_positions.at(vidx);
      
      for (int prev_kidx = 0; prev_kidx < kidx; prev_kidx++) {
        const Vec& kprev_for_vtx = _getkvec(prev_kidx).at(vidx);
        
        Vec kprev_contrib = k_compute_coeffs.at(prev_kidx) * kprev_for_vtx * hstep;
        
        starting_pos = starting_pos + kprev_contrib;
      }
      
      _sys.mesh().vertices().at(vidx).data().r = starting_pos;
    }
    
    // Find which k vector to write into
    vector<double>& kvec = _getkvec(kidx);
    kvec.clear();
    kvec.resize(n_vertices, Vec(0.0, 0.0));
    
    // k vector consists of the velocity evaluated at these starting positions
    // passing by ref to calc instantaneous velocities
    _calc_instantaneous_velocities(kvec, verbose);
  }
  
  void IntegratorAdaptiveRungeKutta::step(bool verbose)
  {
    if (_gamma_param <= 0) {
      throw runtime_error(" IntegratorAdaptiveRungeKutta::step - 'gamma' has not been set to valid value");
    }
    if (_error_allowed_param <= 0) {
      throw runtime_error(" IntegratorAdaptiveRungeKutta::step - 'error_alloweds' has not been set to valid value");
    }
    if (_init_dt <= 0) {
      throw runtime_error(" IntegratorAdaptiveRungeKutta::step - 'init_dt' has not been set to a valid value");
    }
    
    if (_current_dt < 0.0) {
      _current_dt = _init_dt;
    }
    
    _k0.clear();
    _k1.clear();
    _k2.clear();
    _k3.clear();
    _k4.clear();
    _k5.clear();
    _k6.clear();
    
    const double hstep = _current_dt;
    const int n_vertices = _sys.mesh().vertices().size();
    
    // Store the starting positions for reference, before we starting messing around with the vertex data
    const vector<Vec> init_vertex_positions(n_vertices, Vec(0.0, 0.0));
    size_t vidx = 0;
    for (int vidx = 0; vidx < n_vertices; vidx++)
    {
      init_vertex_positions.at(vidx) = _sys.mesh().vertices().at(vidx).data().r;
      // vidx++;
    }
    
    // in order, evaluate each set of k coefficients
    calculate_kvec(0, hstep, init_vertex_positions);
    calculate_kvec(1, hstep, init_vertex_positions);
    calculate_kvec(2, hstep, init_vertex_positions);
    calculate_kvec(3, hstep, init_vertex_positions);
    calculate_kvec(4, hstep, init_vertex_positions);
    calculate_kvec(5, hstep, init_vertex_positions);
    calculate_kvec(6, hstep, init_vertex_positions);
    
    /*
      The fifth order approximation of velocity is the same as k6 - dormand prince is designed this way
    */
    vector<Vec> fifth_order_dr(n_vertices, Vec(0.0,0.0));
    for (int vidx = 0; vidx < n_vertices; vidx++) {
      fifth_order_dr.at(vidx) = _k6.at(vidx) * hstep;
    }
    
    /*
      The fourth order approximation of velocity is calculated using coefficients stored
      in _butcher_tableau_fourth_order_result_coeffs
    */
    vector<Vec> fourth_order_dr(n_vertices, Vec(0.0,0.0));
    for (int vidx = 0; vidx < n_vertices; vidx++) {
      Vec fourth_order_vel = Vec(0.0,0.0);
      
      for (int kidx = 0; kidx < 7; kidx++) {
        double kcoeff = _butcher_tableau_fourth_order_result_coeffs.at(kidx);
        
        fourth_order_vel = fourth_order_vel + ( kcoeff * _getkvec(kidx).at(vidx) );
      }
      
      fourth_order_dr.at(vidx) = fourth_order_vel * hstep;
    }
    
    for (int vidx = 0; vidx < n_vertices; vidx++) {
      double err = fifth_order_dr.at(vidx) - fourth_order_dr.at(vidx);
      
      if (err >= _error_allowed_param) {
        // Reset positions
        for (int vidx = 0; vidx < n_vertices; vidx++) {
          _sys.mesh().vertices().at(vidx).data().r = init_vertex_positions.at(vidx);
        }
        // Set a new, smaller dt
        _current_dt = _calculate_new_timestep(err);
        
        return false;
      }
    }
    
    return true;
  }
}
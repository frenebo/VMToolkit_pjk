/*!
 * \file integrator_adaptive_runge_kutta.cpp
 * \author Paul Kreymborg, frenebo@gmail.com
 * \date 31-Dec-2024
 * \brief IntegratorAdaptiveRungeKutta class 
*/

#include "integrator_adaptive_runge_kutta.hpp"
#include <cmath>


using std::cout;
using std::endl;

namespace VMTutorial
{
  void IntegratorAdaptiveRungeKutta::_calc_instantaneous_velocities(vector<Vec>& vertex_vels_out, bool verbose) const
  {
    double mu = 1.0 / _gamma_param;    // mobility 
    if (verbose || _log_runge_kutta_details) {
      cout << "  calculating instantaneous velocities with mu=1.0/_gamma_param=" << mu << endl;
    }
    
    _force_compute.compute_all_vertex_forces(vertex_vels_out, verbose);
    
    for (auto& vel : vertex_vels_out)
    {
      vel.x *= mu;
      vel.y *= mu;
    }
  }
  
  void IntegratorAdaptiveRungeKutta::_calculate_kvec(int kidx, double hstep, const vector<Vec> & init_positions, bool verbose) {
    // currently unused
    // double dt = hstep * _butcher_tableau_timeoffsets.at(kidx);
    
    const vector<double>& k_compute_coeffs = _butcher_tableau_kcalc_coeffs.at(kidx);
    
    int n_vertices = _sys.mesh().vertices().size();
    
    
    
    // vector<Vec> starting_positions
    for (int vidx = 0; vidx < n_vertices; vidx++) {
      Vec starting_pos = init_positions.at(vidx);
      
      for (int prev_kidx = 0; prev_kidx < kidx; prev_kidx++) {
        Vec kprev_for_vtx = _getkvec(prev_kidx).at(vidx);
        
        Vec kprev_contrib = k_compute_coeffs.at(prev_kidx) * kprev_for_vtx * hstep;
        
        starting_pos = starting_pos + kprev_contrib;
      }
      
      _sys.mesh().vertices().at(vidx).data().r = starting_pos;
    }
    
    // Find which k vector to write into
    vector<Vec>& kvec = _getkvec(kidx);
    kvec.clear();
    kvec.resize(n_vertices, Vec(0.0, 0.0));
    
    // k vector consists of the velocity evaluated at these starting positions
    // passing by ref to calc instantaneous velocities
    _calc_instantaneous_velocities(kvec, verbose);
  }
  
  void IntegratorAdaptiveRungeKutta::adaptive_run_for_time(double t_run,  bool show_progress_bar, bool topological_change, bool verbose)
  {
    const double safety_factor = 0.8; // Multiply the predicted dt by this safety factor so we hopefully undershoot the needed dt
    const int n_vertices = _sys.mesh().vertices().size();
    // If this is the first loop, then we set the current dt to the init dt value;
    if (_current_dt < 0.0) {
      _current_dt = _init_dt;
    }
    
    
    double t_elapsed = 0.0;
    vector<Vec> vtx_dr_out;
    vector<Vec> vtx_errs_out;
    vector<Vec> init_vtx_positions;
    
    while (t_elapsed < t_run) {
      if (topological_change) {
        if (verbose || _log_integrator_stuff) {
          cout << "  IntegratorAdaptiveRungeKutta::adaptive_run_for_time - topological change is enabled, so running topo changes before integrating." << endl;
        }
        
        _topology.T1(verbose);
      }
      double tremaining = t_run - t_elapsed;
      
      
      if (verbose || _log_integrator_stuff){
        cout << "  Storing current vertex positions as starting point for runge-kutta" << endl;
      }
      _get_current_vtx_positions(init_vtx_positions);
      
      bool error_acceptable = false;
      double step_size_final = -1.0;
      
      while (! error_acceptable) {
        if (_current_dt > tremaining || _log_integrator_stuff) {
          if (verbose) {
            cout << "_current_dt is greater than tremaining - " << _current_dt << " > " << tremaining << " : setting current_dt to match tremaining" << endl;
          }
          _current_dt = tremaining;
        }
        _try_step(_current_dt, init_vtx_positions, vtx_dr_out, vtx_errs_out, verbose);
        
        double max_coord_error = 0.0;
        for (int vidx = 0; vidx < n_vertices; vidx++) {
          double verr_x = vtx_errs_out.at(vidx).x;
          double verr_y = vtx_errs_out.at(vidx).y;
          if (verr_x > max_coord_error) max_coord_error = verr_x;
          if (verr_y > max_coord_error) max_coord_error = verr_y;
        }
        
        // Was the error small enough? If so we can move on
        if (max_coord_error > _error_allowed_param) {
          error_acceptable = false;
          
          if (verbose || _log_integrator_stuff) {
            cout << "  Max error (" << max_coord_error << ") was greater than max allowable (" << _error_allowed_param << "), going to retry." << endl;
          }
        } else {
          step_size_final = _current_dt;
          error_acceptable = true;
        }
        
        // After every step, including successful ones, we update the dt stepsize to target the error allowed. this way we aren't stuck with
        // a pointlessly small dt when it is not needed.
        double err_ratio_fifth_root = std::pow( (_error_allowed_param / max_coord_error)  , 1.0/5.0);
        double new_dt = safety_factor * _current_dt * err_ratio_fifth_root;
        if (verbose || _log_integrator_stuff) {
          cout << "  calculating new dt - old dt=" << _current_dt << ", NEW dt=" << new_dt << endl;
        }
        if ( (!error_acceptable) && (new_dt > _current_dt)) {
          throw runtime_error("???????");
        }
        
        _current_dt = new_dt;
      }
      
      // Once the error has been made small enough, we save the vertex displacements and continue
      for (int vidx = 0; vidx < n_vertices; vidx++) {
        _sys.mesh().vertices().at(vidx).data().r += vtx_dr_out.at(vidx);
      }
      
      t_elapsed += step_size_final;
      
      if (show_progress_bar) {
        _update_progress_bar(t_elapsed, t_run, step_size_final);
      }
    }
    
    if (show_progress_bar) {
      _update_progress_bar(t_elapsed, t_run, _current_dt);
    }
    cout << endl;
  }
  
  void IntegratorAdaptiveRungeKutta::_update_progress_bar(double t_elapsed, double t_run, double current_tstep)
  {
    double progress = t_elapsed / t_run;
    
    cout << " [";
    int pos = static_cast<int>(round(_pbar_width * progress));
    for (int i = 0; i < _pbar_width; ++i) {
      if (i < pos) {
        cout << "=";
      } else if (i == pos || i == pos + 1) {
        if (_pbar_spinner_idx == 0) {
          cout << "|";
        } else if (_pbar_spinner_idx == 1) {
          cout << "/";
        } else if (_pbar_spinner_idx == 2) {
          cout << "-";
        } else if (_pbar_spinner_idx == 3) {
          cout << "\\";
        }
      } else {
        cout << " ";
      }
    }
    cout << "] " <<  static_cast<int>(round(progress * 100.0)) << "%" << ", tstep=" << current_tstep << "\r";
    cout.flush();
    
    _pbar_spinner_idx++;
    if (_pbar_spinner_idx > 3) _pbar_spinner_idx = 0;
    
    
  // void Simulation::progress_bar(double progress, const string& end_of_line)
  // {
  //   cout << "[";
  //   int pos = static_cast<int>(round(bar_width * progress));
  //   for (int i = 0; i < bar_width; ++i) 
  //   {
  //     if (i < pos) 
  //       cout << "=";
  //     else if (i == pos) 
  //       cout << ">";
  //     else 
  //       cout << " ";
  //   }
  //   cout << "] "  << static_cast<int>(round(progress * 100.0)) << "%" << end_of_line;
  //   cout.flush();
  // }


    // if (this->print_freq > 0 && !old_style) {
    //   progress_bar(progress, " ");
    // }
    
    // sim_step += steps;
    
    // if (this->print_freq > 0 && !old_style) {
    //   cout << " --> Completed " << sim_step << " simulation steps " << endl;  
    // }
  }
  
  void IntegratorAdaptiveRungeKutta::_try_step(
    double hstep,
    const vector<Vec>& init_vertex_positions,
    vector<Vec>& vtx_dr_out,
    vector<Vec>& vtx_errs_out,
    bool verbose
  ) {
    const int n_vertices = _sys.mesh().vertices().size();
    
    vtx_dr_out.clear(); // We use vtx_dr_out to hold the fourth order dr
    vtx_dr_out.resize(n_vertices, Vec(0.0,0.0));
    vtx_errs_out.clear();
    vtx_errs_out.resize(n_vertices, Vec(0.0,0.0));
    
    /*
      Finds fourth-order and fifth-order runge-kutta results
      returns the fifth-order vertex displacements through vtx_dr_out
      calculates error, returns error through vtx_errs_out
    */
    if (_gamma_param <= 0) {
      throw runtime_error(" IntegratorAdaptiveRungeKutta::step - 'gamma' has not been set to valid value");
    }
    if (_error_allowed_param <= 0) {
      throw runtime_error(" IntegratorAdaptiveRungeKutta::step - 'error_alloweds' has not been set to valid value");
    }
    if (_init_dt <= 0) {
      throw runtime_error(" IntegratorAdaptiveRungeKutta::step - 'init_dt' has not been set to a valid value");
    }
    
    _k0.clear();
    _k1.clear();
    _k2.clear();
    _k3.clear();
    _k4.clear();
    _k5.clear();
    _k6.clear();
    
    // const double hstep = _current_dt;
    
    // in order, evaluate each set of k coefficients
    _calculate_kvec(0, hstep, init_vertex_positions, verbose);
    _calculate_kvec(1, hstep, init_vertex_positions, verbose);
    _calculate_kvec(2, hstep, init_vertex_positions, verbose);
    _calculate_kvec(3, hstep, init_vertex_positions, verbose);
    _calculate_kvec(4, hstep, init_vertex_positions, verbose);
    _calculate_kvec(5, hstep, init_vertex_positions, verbose);
    _calculate_kvec(6, hstep, init_vertex_positions, verbose);
    
    /*
      The fifth order approximation of velocity is the same as k6 - dormand prince is designed this way
    */
    for (int vidx = 0; vidx < n_vertices; vidx++) {
      vtx_dr_out.at(vidx) = _k6.at(vidx) * hstep;
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
      // We use vtx_dr_out to hold the fourth order dr
      // So error is vtx_dr_out - fourth_order_dr
      Vec err = vtx_dr_out.at(vidx) - fourth_order_dr.at(vidx);
      
      vtx_errs_out.at(vidx) = err;
    }
  }
  
  void IntegratorAdaptiveRungeKutta::_get_current_vtx_positions(vector<Vec>& vtx_positions_out) {
    int n_vertices = _sys.mesh().vertices().size();
    
    vtx_positions_out.clear();
    vtx_positions_out.resize(n_vertices, Vec(0.0,0.0));
    for (int vidx=0; vidx < n_vertices; vidx++) {
      vtx_positions_out.at(vidx) = _sys.mesh().vertices().at(vidx).data().r;
    }
  }
  
  
  vector<Vec>& IntegratorAdaptiveRungeKutta::_getkvec(int kidx) {
    if (kidx == 0) return _k0;
    if (kidx == 1) return _k1;
    if (kidx == 2) return _k2;
    if (kidx == 3) return _k3;
    if (kidx == 4) return _k4;
    if (kidx == 5) return _k5;
    if (kidx == 6) return _k6;
    
    throw runtime_error("invalid k idx" + std::to_string(kidx));
  }
}
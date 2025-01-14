/*!
 * \file integrator_adaptive_runge_kutta.hpp
 * \author Paul Kreymborg, frenebo@gmail.com
 * \date 31-Dec-2024
 * \brief IntegratorAdaptiveRungeKutta class 
*/

#ifndef __INTEGRATOR_ADAPTIVE_RUNGE_KUTTA_HPP__
#define __INTEGRATOR_ADAPTIVE_RUNGE_KUTTA_HPP__


#include "integrator.hpp"
#include <map>
#include <stdexcept>

using std::map;
using std::runtime_error;
using std::cout;
using std::endl;


namespace VMTutorial
{

  class IntegratorAdaptiveRungeKutta : public Integrator 
  {
    /*
      Implementation of Dormand-Prince Runge Kutta, with adaptive timestep.
    */
    public:
      
      IntegratorAdaptiveRungeKutta(
        System& sys,
        ForceCompute& fc,
        Topology& topology
      ) : Integrator{
            sys,
            fc,
            topology,
            Integrator::IntegrationType::ADAPTIVE_TIMESTEP
          },
          _gamma_param{-1.0},
          _error_allowed_param{-1.0},
          _init_dt{-1.0},
          _current_dt{-1.0},
          _log_integrator_stuff{false},
          _log_runge_kutta_details{false},
          _pbar_spinner_idx{0},
          _pbar_width{40},
          
          
          // Butcher talbeau for the Runge-Kutta Dormand-Prince method
          _butcher_tableau_timeoffsets{
            0,        // k0
            1.0/5.0,  // k1
            3.0/10.0, // k2
            4.0/5.0,  // k3
            8.0/9.0,  // k4
            1.0,      // k5
            1.0,      // k6
          },
          _butcher_tableau_kcalc_coeffs{
            {}, // k0 is f(t_init, y_init + 0)
            {1.0/5.0}, // k1 is f(t_init + 0.2*tstep, y_init + 0.2*k0*tstep)
            {3.0/40.0,       9.0/40.0}, // k2 coeffs
            {44.0/45.0,      -56.0/15.0,      32.0/9.0}, // k3 coeffs
            {19372.0/6561.0, -25360.0/2187.0, 64448.0/6561.0, -212.0/729.0}, // k4 coeffs
            {9017.0/3168.0,  -355.0/33.0,     46732.0/5247.0, 49.0/176.0,  -5103.0/18656.0}, // k5 coeffs
            {35.0/384.0,     0.0,             500.0/1113.0,   125.0/192.0, -2187.0/6784.0,    11.0/84.0} // k6 coeffs
          },
          // fifth order is calculated directly from the k coeffs - coefficients are identical to the k6 coeffs, so no need to recalculate
          // fourth order approx is y_init + (some val)*k0*tstep + etc.
          _butcher_tableau_fourth_order_result_coeffs{
            5179.0/57600.0,  0.0,             7571.0/16695.0, 393.0/640.0, -92097.0/339200.0, 187.0/2100.0, 1.0/40.0
          }
      {
        cout << "IntegratorAdaptiveRungeKutta - NOTE: this implementation of dormand-prince RK has NOT been verified thoroughly - the error calculated may not be correct." << endl;
        // blah blah
      }
      
      void adaptive_run_for_time(double t_run, bool show_progress_bar, bool topological_change, bool verbose) override;
      
      void set_params(const params_type& params) override 
      { 
        for (auto& p : params)
        {
          if (p.first == "gamma")
          {
            _gamma_param = p.second;
          } else if (p.first == "error_allowed") {
            _error_allowed_param = p.second;
          } else if (p.first == "init_dt") {
            _init_dt = p.second;
          } else {
            throw runtime_error("Unknown parameter name - " + p.first);
          }
        }
      }
      
      void set_flag(const string& flagname, bool flag_val) override
      {
        if (flagname == "log_integrator_details") {
          _log_integrator_stuff = flag_val;
        } else if (flagname == "log_runge_kutta_details") {
          _log_runge_kutta_details = flag_val;
        } else {
          throw runtime_error("IntegratorAdaptiveRungeKutta::set_flag - unknown flagname " + flagname);
        }
      }
      
      
    private:
      void _calc_instantaneous_velocities(vector<Vec>& vertex_vels_out, bool verbose) const;
      vector<Vec>& _getkvec(int kidx);
      void _calculate_kvec(int kidx, double hstep, const vector<Vec> & init_positions, bool verbose);
      void _get_current_vtx_positions(vector<Vec>& vtx_positions_out);
      void _try_step(
        double hstep, // time step size
        const vector<Vec>& init_vertex_positions,
        vector<Vec>& vtx_dr_out,
        vector<Vec>& vtx_errs_out,
        bool verbose
      );
      
      void _update_progress_bar(double t_elapsed, double t_run, double current_tstep);
      
      double _gamma_param;
      double _error_allowed_param;
      double _init_dt;
      double _current_dt;
      bool _log_integrator_stuff;
      bool _log_runge_kutta_details;
      int _pbar_spinner_idx;
      int _pbar_width;
      
      // These k values are cleared and recalculated every timestep as part of the runge-kutta algorithm
      vector<Vec> _k0;
      vector<Vec> _k1;
      vector<Vec> _k2;
      vector<Vec> _k3;
      vector<Vec> _k4;
      vector<Vec> _k5;
      vector<Vec> _k6;
      
      const vector<double> _butcher_tableau_timeoffsets;
      const vector<vector<double>> _butcher_tableau_kcalc_coeffs;
      const vector<double> _butcher_tableau_fourth_order_result_coeffs;
  };
}


#endif

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


namespace VMTutorial
{

  class IntegratorAdaptiveRungeKutta : public Integrator 
  {
    /*
      Implementation of RK4
    */
    public:
      
      IntegratorAdaptiveRungeKutta(
        System& sys,
        ForceCompute& fc
      ) : Integrator{sys, fc},
          _gamma_param{-1.0},
          _error_allowed_param{-1.0},
          _init_dt{-1.0},
          _current_dt{-1.0}
      {
        // Runge-Kutta Dormand-Prince method
        _butcher_tableau_timeoffsets = vector<double> {
          0,        // k0
          1.0/5.0,  // k1
          3.0/10.0, // k2
          4.0/5.0,  // k3
          8.0/9.0,  // k4
          1.0,      // k5
          1.0,      // k6
        };
        
        _butcher_tableau_kcalc_coeffs = vector<vector<double>>  {
          {}, // k0 is f(t_init, y_init + 0)
          {1.0/5.0}, // k1 is f(t_init + 0.2*tstep, y_init + 0.2*k0*tstep)
          {3.0/40.0,       9.0/40.0}, // k2 coeffs
          {44.0/45.0,      -56.0/15.0,      32.0/9.0}, // k3 coeffs
          {19372.0/6561.0, -25360.0/2187.0, 64448.0/6561.0, -212.0/729.0}, // k4 coeffs
          {9017.0/3168.0,  -355.0/33.0,     46732.0/5247.0, 49.0/176.0,  -5103.0/18656.0}, // k5 coeffs
          {35.0/384.0,     0.0,             500.0/1113.0,   125.0/192.0, -2187.0/6784.0,    11.0/84.0} // k6 coeffs
        };
        // fifth order is calculated directly from k coeffs - coefficients are identical
        // fourth order approx is y_init + (some val)*k0*tstep + etc.
        _butcher_tableau_fourth_order_result_coeffs = vector<double> {
          5179.0/57600.0,  0.0,             7571.0/16695.0, 393.0/640.0, -92097.0/339200.0, 187.0/2100.0, 1.0/40.0
        };
      }
      
      void step(bool verbose) override;
      
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
      
      
    private:
      vector<Vec> _instantaneous_velocities(bool verbose) const;
      
      vector<double> _k0;
      vector<double> _k1;
      vector<double> _k2;
      vector<double> _k3;
      vector<double> _k4;
      vector<double> _k5;
      vector<double> _k6;
      
      vector<double> _butcher_tableau_timeoffsets;
      vector<vector<double>> _butcher_tableau_kcalc_coeffs;
      vector<double> _butcher_tableau_fourth_order_result_coeffs;
      // vector<double> _butcher_tableau_fifth_order_result_coeffs;
      double _gamma_param;
      double _error_allowed_param;
      double _init_dt;
      double _current_dt;
  };
}


#endif

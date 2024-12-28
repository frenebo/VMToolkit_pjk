/*!
 * \file integrator_runge_kutta.cpp
 * \author Paul Kreymborg, frenebo@gmail.com
 * \date 27-Dec-2024
 * \brief IntegratorRungeKutta class 
*/

#include  "integrator_runge_kutta.hpp"


namespace VMTutorial
{
  void IntegratorRungeKutta::step(bool verbose)
  {
    if (verbose)
    {
        cout << "IntegratorRungeKutta::step - executing with dt=" << _dt << endl;
    }
    
    
    
    throw runtime_error("IntegratorRungeKutta::step -  Unimplemented");
  }
}
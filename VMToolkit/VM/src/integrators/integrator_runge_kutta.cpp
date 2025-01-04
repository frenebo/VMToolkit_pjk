/*!
 * \file integrator_runge_kutta.cpp
 * \author Paul Kreymborg, frenebo@gmail.com
 * \date 27-Dec-2024
 * \brief IntegratorRungeKutta class 
*/

#include <vector>

#include  "integrator_runge_kutta.hpp"

using std::cout;
using std::vector;
using std::endl;

namespace VMTutorial
{
  vector<Vec> IntegratorRungeKutta::_instantaneous_velocities(bool verbose) const
  {
    double mu = 1.0 / _gamma;    // mobility 
    if (verbose) {
      cout << "  calculating instantaneous velocities with mu=1.0/_gamma=" << mu << endl;
    }
    
    // verbose = false;
    
    vector<Vec> vertex_vels;
    _force_compute.compute_all_vertex_forces(vertex_vels, verbose);
    
    // v = F/gamma
    for (auto& vel : vertex_vels)
    {
      vel.x *= mu;
      vel.y *= mu;
    }
    
    return vertex_vels;
  }
  
  void IntegratorRungeKutta::step(bool verbose)
  {
    if (verbose)
    {
        cout << "IntegratorRungeKutta::step - executing with dt=" << _dt << endl;
    }
    
    double h = _dt;
    
    size_t n_vertices = _sys.mesh().vertices().size();
    vector<Vec> init_vertex_positions(n_vertices, Vec(0.0, 0.0));
    
    size_t vidx = 0;
    for (const auto& v : _sys.mesh().vertices())
    {
      init_vertex_positions.at(vidx) = v.data().r;
      vidx++;
    }
    
    vector<Vec> k1 = instantaneous_velocities(verbose); // Vertex velocities starting from init conditionss
    
    // k2 = f(y_n + h*k1*1/2)
    vidx = 0;
    for (auto& v : _sys.mesh().vertices()) {
      Vec dr_fork2 = h * k1.at(vidx) * 0.5;
      v.data().r = init_vertex_positions.at(vidx) + dr_fork2;
      
      if (verbose && vidx==0){ cout << "    dr_fork2[0]="<<dr_fork2<<endl; }
      
      vidx++;
    }
    vector<Vec> k2 = instantaneous_velocities(verbose);
    
    // k3 = f(y_n + h*k2 * 1/2)
    vidx = 0;
    for (auto& v : _sys.mesh().vertices()) {
      Vec dr_fork3 = h * k2.at(vidx) * 0.5;
      v.data().r = init_vertex_positions.at(vidx) + dr_fork3;
      
      if (verbose && vidx==0){ cout << "    dr_fork3[0]="<<dr_fork3<<endl; }
      
      vidx++;
    }
    vector<Vec> k3 = instantaneous_velocities(verbose);
    
    // k4 = f(y_n + h*k3)
    vidx = 0;
    for (auto& v : _sys.mesh().vertices()) {
      Vec dr_fork4 = h * k3.at(vidx);
      v.data().r = init_vertex_positions.at(vidx) + dr_fork4;
      
      
      if (verbose && vidx==0){ cout << "    dr_fork4[0]="<<dr_fork4<<endl; }
      
      vidx++;
    }
    vector<Vec> k4 = instantaneous_velocities(verbose);
    
    // y_(n+1) = y_n + (h/6)*( k1 + 2*k2 + 2*k3 + k4 )
    vidx = 0;
    double dt_inv = 1/_dt;
    for (auto& v : _sys.mesh().vertices()) {
      Vec dr = (h/6) * (
        k1.at(vidx) + 2.0*k2.at(vidx) + 2.0*k3.at(vidx) + k4.at(vidx)
      );
      
      Vec vel_new = dr * dt_inv;
      
      v.data().r = init_vertex_positions.at(vidx) + dr;
      v.data().vel = vel_new;
      
      vidx++;
    }
    
    if (verbose)
    {
      cout << " Finished step - k1[0],k2[0],k3[0],k4[0] : " << k1.at(0) << k2.at(0) << k3.at(0) << k4.at(0) << endl;
    }
  }
}

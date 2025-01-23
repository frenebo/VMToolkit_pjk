/*!
 * \file constraint_manager.hpp
 * \author Paul Kreymborg, pk5192@princeton.edu
 * \date 09-Jan-2025
 * \brief ConstraintManager class
 */
 
#ifndef __CONSTRAINT_MANAGER_HPP__
#define __CONSTRAINT_MANAGER_HPP__

#include <vector>
#include <string>

#include "../system.hpp"
#include <stdexcept>


namespace VMSim
{
  using std::runtime_error;

  class ConstraintManager: public ClassFactory<Constrainer>
  {
    ConstraintManager(System& sys) : _sys{sys}
    {
    }
    
    ~ConstraintManager() = default; 
      
    private:
      System& _sys;
  }
}


#endif

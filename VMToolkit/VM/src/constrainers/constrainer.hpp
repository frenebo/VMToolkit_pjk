/*!
 * \file constrainer.hpp
 * \author Paul Kreymborg, pk5192@princeton.edu
 * \date 09-Jan-2025
 * \brief Constrainer class
 */
 
#ifndef __CONSTRAINER_HPP__
#define __CONSTRAINER_HPP__

#include "../system.hpp"


namespace VMTutorial
{
  class Constrainer
  {
    
      Integrator(System& sys) : _sys{sys}
      { 
        
      }
      virtual ~Integrator() { }
  }
}
#endif
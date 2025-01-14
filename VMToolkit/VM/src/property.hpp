/*!
 * \file property.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 31-Nov-2023
 * \brief Property class
 */

#ifndef __PROPERTY_HPP__
#define __PROPERTY_HPP__

#include <map>
#include <vector>
#include <string>

#include "vec.hpp"


namespace VMSim
{

	// struct Property
	// {
	struct HEProperty
	{
	};
	struct VertexProperty
	{
		Vec vel;
		Vec r;
		
		// VertexProperty &operator=(const VertexProperty &p)
		// {
		// 	if (this == &p)
		// 		return *this;
			
		// 	return *this;
		// }
	};
	struct EdgeProperty
	{
	};
	struct FaceProperty
	{
		// FaceProperty &operator=(const FaceProperty &p)
		// {
		// 	if (this == &p)
		// 		return *this;
			
		// 	return *this;
		// }
	};
	// };

}

#endif

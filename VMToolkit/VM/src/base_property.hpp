
/*!
 * \file base_property.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 30-Nov-2023
 * \brief BaseProperty class
 */

#ifndef __BASE_PROPERTY_HPP__
#define __BASE_PROPERTY_HPP__

namespace VMTutorial
{
	struct BaseProperty
	{
		struct HEProperty
		{
		};
		struct VertexProperty
		{
			VertexProperty &operator=(const VertexProperty &p)
			{
				if (this == &p)
					return *this;
				
				return *this;
			}
		};
		struct EdgeProperty
		{
		};
		struct FaceProperty
		{
			FaceProperty &operator=(const FaceProperty &p)
			{
				if (this == &p)
					return *this;
				
				return *this;
			}
		};
	};
}

#endif

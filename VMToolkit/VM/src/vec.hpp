/*!
 * \file vec.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 30-Nov-2023
 * \brief 2D vectors
 */

#ifndef __VEC_HPP__
#define __VEC_HPP__

#include <cmath>
#include <memory>
#include <iostream>

#include <iomanip>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>


namespace py = pybind11;

namespace VMSim
{

	class Vec
	{
	public:
		//! Default constructor
		Vec() : x{0.0}, 
				y{0.0} 
		{

		}
		//! Constructor for a vector object
		Vec(double x, double y) : x{x}, 
								  y{y}  
		{

		}
		
		//! Copy constructor
		Vec(const Vec &v)
		{
			x = v.x;
			y = v.y;
		}
		//! Assignment operator
		Vec &operator=(const Vec &rhs)
		{
			x = rhs.x;
			y = rhs.y;
			
			return *this;
		}

		//! Add two vectors
		Vec operator+(const Vec &v)
		{
			double xx = x + v.x, yy = y + v.y;
			
			return Vec(xx, yy);
		}

		//! Add two vectors (constant version)
		Vec operator+(const Vec &v) const { return *(const_cast<Vec *>(this)) + v; }

		//! Subtract two vectors
		Vec operator-(const Vec &v)
		{
			double xx = x - v.x, yy = y - v.y;
			
			return Vec(xx, yy);
		}

		//! Subtract two vectors (constant version)
		Vec operator-(const Vec &v) const { return *(const_cast<Vec *>(this)) - v; }

		//! Negate a vector
		Vec operator-()
		{
			return Vec(-x, -y);
		}

		//! Scale vector by a constant
		Vec operator*(const double c)
		{
			double xx = c * x, yy = c * y;
			
			return Vec(xx, yy);
		}

		//! Scale vector by a constant (const version)
		Vec operator*(const double c) const { return (*(const_cast<Vec *>(this))) * c; }

		//! Test equality
		bool operator==(const Vec &v)
		{
			return (x == v.x && y == v.y);
		}

		//! Add vector to current vector
		Vec &operator+=(const Vec &v)
		{
			double xx = x + v.x, yy = y + v.y;
			
			x = xx;
			y = yy;
			return *this;
		}

		//! Subtract vector from current vector
		Vec &operator-=(const Vec &v)
		{
			double xx = x - v.x, yy = y - v.y;
			
			x = xx;
			y = yy;
			return *this;
		}

		//! Euclidean dot product with another vector
		double dot(const Vec &v)
		{
			return x * v.x + y * v.y;
		}

		//! Vector length
		double len() const { return std::sqrt(x * x + y * y); }

		//! Vector length squared
		double len2() { return x * x + y * y; }

		//! Make the vector has unit length
		void normalise()
		{
			double len = this->len();
			if (len != double(0))
			{
				x /= len;
				y /= len;
			}
		}

		//! Return unit vector in the direction of this vector
		Vec unit()
		{
			double len = this->len();
			if (len != double(0))
				return Vec(x / len, y / len);
			return Vec(x, y);
		}

		Vec unit() const
		{
			double len = this->len();
			if (len != double(0))
				return Vec(x / len, y / len);
			return Vec(x, y);
		}

		//! Rotate vector by \f$ \phi \f$
		Vec rotate(const double phi)
		{
			double s = std::sin(phi);
			double c = std::cos(phi);

			double xx = c * x - s * y;
			double yy = s * x + c * y;


			return Vec(xx, yy);
		}

		//! Compute e_z x v (used for force)
		Vec ez_cross_v() { return Vec(-y, x); }

		Vec ez_cross_v() const { return Vec(-y, x); }


		// Scale vector
		void scale(double a, double b)
		{
			this->x = a * x;
			this->y = b * y;
		}

		// friend Vec operator*(const Matrix &, const Vec &);

		double x, y; //!< Position in the embedding 3d flat space

	};

	//! Scale vector by a number
	inline Vec operator*(const double c, const Vec &v)
	{
		double xx = c * v.x, yy = c * v.y;
		
		return Vec(xx, yy);
	}

	// inline Vec operator*(const Matrix &m, const Vec &v)
	// {
	// 	return Vec(m._mxx * v.x + m._mxy * v.y, m._myx * v.x + m._myy * v.y);
	// }

	// inline double dot(const Vec &v1, const Vec &v2)
	// {
	// 	return (v1.x * v2.x + v1.y * v2.y);
	// }

	// inline double cross(const Vec &v1, const Vec &v2)
	// {
	// 	return (v1.x * v2.y) - (v1.y * v2.x);
	// }

	inline std::ostream &operator<<(std::ostream &os, const Vec &v)
	{
		os << std::setprecision(22) << "(" << v.x << "," << v.y << ")";
		return os;
	}

	void export_Vec(py::module &m);

}

#endif

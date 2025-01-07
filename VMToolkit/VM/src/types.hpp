/*!
 * \file types.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 30-Nov-2023
 * \brief Data types for the Mesh
 */

#ifndef __TYPES_HPP__
#define __TYPES_HPP__

// #define VERSION

#include <vector>
#include <iterator>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "vec.hpp"
#include "property.hpp"

namespace py = pybind11;

namespace VMTutorial
{
	// Circulator definitions
	class VertexCirculator;
	class FaceCirculator;
	class VertexCCirculator;
	class FaceCCirculator;
	
	// Forward declarations
	class HalfEdge;
	
	class Vertex;
	
	class Edge;
	
	class Face;
	
	class Mesh;
	// Handles to mesh elements
	
	using HEHandle = typename std::vector<HalfEdge>::iterator;
	
	using VertexHandle = typename std::vector<Vertex>::iterator;
	
	using EdgeHandle = typename std::vector<Edge>::iterator;
	
	using FaceHandle = typename std::vector<Face>::iterator;

	
	using HECHandle = typename std::vector<HalfEdge>::const_iterator;
	
	using VertexCHandle = typename std::vector<Vertex>::const_iterator;
	
	using EdgeCHandle = typename std::vector<Edge>::const_iterator;
	
	using FaceCHandle = typename std::vector<Face>::const_iterator;
	

	//! HalfEdge class
	class HalfEdge
	{

	public:
		// Constructors
		HalfEdge(Mesh &mesh) : _mesh{mesh},
										 _idx{-1},
										 _property{HEProperty()},
										 _from{-1},
										 _to{-1},
										 _edge{-1},
										 _face{-1},
										 _pair{-1},
										 _next{-1},
										 _prev{-1}
										//  erased{false}
		{
		}

		// Member functions
		int idx() const { return _idx; }
		void set_idx(int idx) { _idx = idx; }

		HEProperty &data() { return _property; }
		HEProperty data() const { return _property; }

		VertexHandle from();
		VertexHandle to();

		VertexCHandle from() const;
		VertexCHandle to() const;

		EdgeHandle edge();
		EdgeCHandle edge() const;

		FaceHandle face();
		FaceCHandle face() const;

		HEHandle pair();
		HEHandle next();
		HEHandle prev();

		HECHandle pair() const;
		HECHandle next() const;
		HECHandle prev() const;
		
		void set_pair(int he_idx) { _pair = he_idx; }
		void set_next(int he_idx) { _next = he_idx; }
		void set_prev(int he_idx) { _prev = he_idx; }
		
		void set_from(int vtx_idx) { _from = vtx_idx; }
		void set_to(int vtx_idx) { _to = vtx_idx; }
		
		void set_face(int face_idx) { _face = face_idx; }

		Vec direction();

		friend class Mesh; // This is used to we can read in the mesh without too much boilerplate code

	private:
		int _idx; // Half-edge index (used for debugging)
		HEProperty _property;

		int _from; // vertex it starts from (that vertex will have this halfedge as its he)
		int _to;   // vertex it points to (its pair will have this vertex as its he)

		int _edge; // edge this he is part of

		int _face; // face to the left of it, when looking in the direction this he points to

		int _pair; // its pair half edge (undefined for boundary edges)
		int _next; // next he in the same face
		int _prev; // previous he in the same face

		Mesh &_mesh; // This can be avoided if one uses lists and list iterators
	};

	//!< Vertex class
	class Vertex
	{
	public:
		// Constructors
		Vertex(Mesh &mesh) : _mesh{mesh},
									   id{0},
									   _he{-1},
									   _property{VertexProperty()},
									   boundary{false}
		{
			_property.r = Vec(0.0,0.0);
		}
		Vertex(int id, const Vec &r, Mesh &mesh) : _mesh{mesh},
															 id{id},
															 _property{VertexProperty()},
															 _he{-1},
															//  erased{false},
															 boundary{false}
		{
			_property.r = r;
		}
		Vertex(int id, const Vec &r, bool bnd, Mesh &mesh) : _mesh{mesh},
																	   id{id},
																	   _he{-1},
																	   _property{VertexProperty()},
																	//    erased{false},
																	   boundary{bnd}
		{
			_property.r = r;
		}

		
		VertexProperty &data();
		VertexProperty data() const;

		HEHandle he();
		HECHandle he() const;
		
		void set_he(int he_idx) { _he = he_idx; }
		
		// void idx() { return id;}

		VertexCirculator circulator();
		VertexCCirculator circulator() const;

		// Public members
		// Vec r;			  // position
		int id;			  // unique id
		// bool erased;	  // marks vertices that are not connected to the rest of the mesh, but are still in memory
		bool boundary;	  // if true, vertex is on boundary
		int coordination; // number of neighbours this vertex has

		friend class Mesh; // This is used to we can read in the mesh without too much boilerplate code

	private:
		VertexProperty _property;
		int _he; // outgoing half edge
		Mesh &_mesh;
	};

	//!< Edge class
	class Edge
	{
	public:
		// Constructors
		Edge(Mesh &mesh) : _mesh{mesh},
									 _idx{0},
									 i{0},
									 j{0},
									 _he{-1},
									 _property{EdgeProperty()},
									 boundary{true}
		{
		}
		Edge(int i, int j, Mesh &mesh) : _mesh{mesh},
												   _idx{0},
												   i{i},
												   j{j},
												   _property{EdgeProperty()},
												   boundary{true}
		{
		}

		// Member functions
		int idx() const { return _idx; }
		void set_idx(int idx) { _idx = idx; }
		EdgeProperty &data() { return _property; }
		EdgeProperty data() const { return _property; }

		HEHandle he();
		HECHandle he() const;

		// Public members
		int i, j;	   // indices of two vertices
		bool boundary; // if true, edge is a boundary edge
		// bool erased;   // marks all erased edges that are still in memory

		friend class Mesh;

	private:
		int _idx; // Unique edge index
		EdgeProperty _property;
		int _he; // one of the two half edges
		Mesh &_mesh;
	};

	//!< Face class
	class Face
	{
	public:
		// Constructors
		Face(Mesh &mesh) : _mesh{mesh},
									 id{0},
									 _he{-1},
									 _property{FaceProperty()},
									 outer{false}
		{
		}
		Face(int id, Mesh &mesh) : _mesh{mesh},
											 id{id},
											 _he{-1},
											 _property{FaceProperty()},
											 outer{false}
		{
		}

		// Member functions
		FaceProperty &data() { return _property; }
		FaceProperty data() const { return _property; }

		HEHandle he();
		HECHandle he() const;
		
		void set_he(int he_idx) { _he = he_idx; }

		FaceCirculator circulator();
		FaceCCirculator circulator() const;

		// public members
		int id;		 // face id
		bool outer;	 // if true, face is a ghost outer face
		int nsides;	 // number of sides face has
		// bool erased; // if true, face is marked as erased

		friend class Mesh;

	private:
		FaceProperty _property;
		int _he; // one of its half edges
		Mesh &_mesh;
	};
}

#endif

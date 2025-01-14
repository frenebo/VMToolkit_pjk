
/*!
 * \file mesh.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 09-Jun-2017
 * \brief Mesh class
 */

#ifndef __MESH_HPP__
#define __MESH_HPP__

#include <vector>
#include <iostream>
#include <string>

#include "types.hpp"
#include "type_circulators.hpp"

using std::runtime_error;
using std::vector;
using std::string;

namespace VMSim
{

	class Mesh
	{
	public:
		Mesh() : _log_topology_stuff{false}
		{
			//
		}

		//! Mesh setup functions
		void add_vertex(const Vertex &v)
		{
			_vertices.push_back(v);
		}
		void add_edge(const Edge &e) { _edges.push_back(e); }
		void add_halfedge(const HalfEdge &he) { _halfedges.push_back(he); }
		void add_face(int face_id, const vector<int> &vert_ids, bool verbose=false);


		void wipe()
		{
			_halfedges.clear();
			_vertices.clear();
			_edges.clear();
			_faces.clear();
		}

		// accessor functions
		HEHandle get_he(int);
		VertexHandle get_vertex(int);
		EdgeHandle get_edge(int);
		FaceHandle get_face(int);

		Vertex &get_vertex_reference(int);
		// HalfEdge &get_halfedge(int);
		// Edge &get_edge(int);
		Face &get_face_reference(int);

		vector<HalfEdge> &halfedges() { return _halfedges; }
		vector<Vertex> &vertices() { return _vertices; }
		const vector<Vertex> &cvertices() const { return _vertices; }
		vector<Edge> &edges() { return _edges; }
		const vector<Edge> &cedges() const { return _edges; }
		vector<Face> &faces() { return _faces; }
		const vector<Face> &cfaces() const { return _faces; }

		int num_vert() const { return _vertices.size(); }
		int num_faces() const { return _faces.size(); }

		// mesh manipulation functions
		bool T1(Edge &, double, bool verbose);

		// mesh info functions
		// @TODO move to Dumps class
		vector<vector<double>> get_vertex_positions() const
		{
			vector<vector<double>> vpositions(num_vert(), vector<double>{0.0,0.0});
			
			size_t vidx = 0;
			for (const auto& v : _vertices)
			{
				vpositions[vidx][0] = v.data().r.x;
				vpositions[vidx][1] = v.data().r.y;
				
				vidx++;
			}
			
			return vpositions;
		}
		
		// @TODO move to Dumps class
		vector<vector<int>> get_face_member_vertex_ids() const
		{
			vector<vector<int>> faces_members;
			
			size_t fidx = 0;
			for (const auto& f : _faces)
			{
				vector<int> vids;
				for (auto he : f.circulator())
				{
					vids.push_back(he.from()->id);
				}
				faces_members.push_back(vids);
			}
			
			return faces_members;
		}

		double area(const Face &) const;
		double perim(const Face &) const;
		double len(const Edge &);
		int coordination(const Vertex &);
		int face_sides(const Face &);
		bool is_boundary_face(const Face &);

		void tidyup();	  // get the mesh in order (set boundary edges, outer faces, etc. )
		// Vec get_centre(); // compute geomes
		Vec get_face_centroid(const Face &);
		
		void set_logflag(string logname, bool new_val) {
			if (logname == "topology") {
				_log_topology_stuff = new_val;
			} else {
				throw runtime_error("Unknown log flag name " + logname);
			}
		}

	private:
		bool _edge_allowed_to_transition(const Edge &e, bool verbose) const;
		
		bool _log_topology_stuff;
		
		std::vector<HalfEdge> _halfedges;
		std::vector<Vertex> _vertices;
		std::vector<Edge> _edges;
		std::vector<Face> _faces;
		
	};

};

#endif

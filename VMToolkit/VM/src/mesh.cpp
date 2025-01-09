
#include "mesh.hpp"

using std::runtime_error;
using std::vector;
using std::endl;
using std::cout;

namespace VMTutorial
{
    
	Vertex &Mesh::get_vertex_reference(int i)
	{
		if ((i < 0) || (i > _vertices.size()))
			throw runtime_error("Vertex index out of bounds.");
		else
			return _vertices[i];
	}

	// HalfEdge &Mesh::get_halfedge(int i)
	// {
	// 	if ((i < 0) || (i > _halfedges.size()))
	// 		throw runtime_error("Junction index out of bounds.");
	// 	else
	// 		return _halfedges[i];
	// }

	// Edge &Mesh::get_edge(int i)
	// {
	// 	if ((i < 0) || (i > _edges.size()))
	// 		throw runtime_error("Junction index out of bounds.");
	// 	else
	// 		return _edges[i];
	// }

	Face &Mesh::get_face_reference(int i)
	{
		if ((i < 0) || (i > _faces.size()))
			throw runtime_error("Face index out of bounds.");
		else
			return *(std::find_if(_faces.begin(), _faces.end(), [i](const Face &f) -> bool
							 { return (f.id == i); }));
	}

	// accessor functions
	HEHandle Mesh::get_he(int i)
	{
		if ((i >= 0) && (i < _halfedges.size()))
		{
			return std::next(_halfedges.begin(), i);
		}
		else
		{
			throw runtime_error("HalfEdge index out of bounds.");
		}
	}

	VertexHandle Mesh::get_vertex(int i)
	{
		if ((i >= 0) && (i < _vertices.size()))
		{
			return std::next(_vertices.begin(), i);
		}
		else
		{
			throw runtime_error("Vertex index out of bounds.");
		}
	}

	EdgeHandle Mesh::get_edge(int i)
	{
		if ((i >= 0) && (i < _edges.size()))
		{
			return std::next(_edges.begin(), i);
		}
		else
		{
			throw runtime_error("Edge index out of bounds.");
		}
	}

	FaceHandle Mesh::get_face(int i)
	{
		if ((i >= 0) && (i < _faces.size()))
		{
			return std::next(_faces.begin(), i);
		}
		else
		{
			throw runtime_error("Face index out of bounds.");
		}
	}
	// end of accessor functions

	// mesh setup functions
	void Mesh::add_face(int face_id, const vector<int> &vert_ids, bool verbose)
	{
		if (verbose) { std::cout << "Adding face" << std::endl; }
		_faces.push_back(Face(face_id, *this));
		
		FaceHandle fh = std::prev(_faces.end());
		HEHandle he;
		int prev_he;
		int first_he;
		if (verbose) {
			std::cout << "Adding vertices" << std::endl;
			std::cout << "  # of vertice in vert_ids: " << vert_ids.size() << std::endl;
		}
		for (int i = 0; i < vert_ids.size(); i++)
		{
			if (verbose) { std::cout << "  vtx idx " << i << std::endl; }
			int v1_id = vert_ids[i], v2_id = vert_ids[(i == (vert_ids.size() - 1)) ? 0 : i + 1];
			
			VertexHandle vh_from = this->get_vertex(v1_id);
			
			VertexHandle vh_to = this->get_vertex(v2_id);
			
			
			HalfEdge he_new = HalfEdge(*this);
			
			_halfedges.push_back(he_new);
			
			_halfedges.back().set_idx(_halfedges.size() - 1);
			
			he = this->get_he(_halfedges.size() - 1);
			
			if (i == 0) {
				first_he = he->idx();
			}
			he->_from = vh_from->id;
			he->_to = vh_to->id;
			vh_from->_he = he->idx();
			EdgeHandle eh = std::find_if(_edges.begin(), _edges.end(), [v1_id, v2_id](const Edge &e) -> bool
												{ return (v1_id == e.i && v2_id == e.j) || (v1_id == e.j && v2_id == e.i); });
			
			if (eh == _edges.end())
			{
				_edges.push_back(Edge(v1_id, v2_id, *this));
				_edges.back().set_idx(_edges.size() - 1);
				eh = std::prev(_edges.end());
				eh->_he = he->idx();
			}
			else
			{
				eh->he()->_pair = he->idx();
				he->_pair = eh->he()->idx();
				eh->boundary = false;
			}
			he->_edge = eh->idx();
			he->_face = fh->id;
			
			if (i > 0)
			{
				he->_prev = _halfedges[prev_he].idx();
				_halfedges[prev_he]._next = he->idx();
			}
			prev_he = he->idx();
		}
		_halfedges[first_he]._prev = he->idx();
		he->_next = _halfedges[first_he].idx();
		fh->_he = _halfedges[first_he].idx(); // Make sure that the first half edge is the face half-edge. This makes life easier when reading in "per edge" data.
		fh->nsides = this->face_sides(*fh);
	}
	
	bool Mesh::_edge_allowed_to_transition(const Edge &e, bool verbose) const
	{
		// In order to transition, the edge must have different faces on each side (duh),
		// and also two different faces that it's vertices intersect - so that
		// after the t1 transition, it will have different faces on each side.
		// This excludes edges whose neighbor has the exact same neighboring faces.
		
		HECHandle he = e.he();		 // Half-edge belonging to the end
		HECHandle hep = he->pair(); // Half-edge pair to he
		
		
		// There must be a 'branch' after the edge
		int idx_he_next_pair = he->next()->pair()->idx();
		int idx_hep_prev = hep->prev()->idx();
		
		if (idx_he_next_pair == idx_hep_prev) {
			if (verbose) {
				cout << "  edge's next edge does not branch - not allowed to transition" << endl;
			}
			return false;
		}
		
		// There must be a 'branch' before the edge
		int idx_he_prev_pair = he->prev()->pair()->idx();
		int idx_hep_next = hep->next()->idx();
		
		if (idx_he_prev_pair == idx_hep_next) {
			if (verbose) {
				cout << "  edge's prev edge does not branch - not allowed to transition" << endl;
			}
			return false;
		}
		
		return true;
	}

	// Mesh manipulation functions
	//! Implements actual T1
	bool Mesh::T1(Edge &e, double edge_len, bool verbose)
	{
		if (! _edge_allowed_to_transition(e, verbose)) {
			if (verbose || _log_topology_stuff) {
				cout << "    Mesh::T1 - edge not allowed to transition, skipping T1 transition" << endl;
			}
			return false;
		}
		

		// he->face()->nsides = this->face_sides(*(he->face()));
		// hep->face()->nsides = this->face_sides(*(hep->face()));
		
		// cout << "finished setting nsides" << endl;
		

		HEHandle he = e.he();		 // Half-edge belonging to the end
		HEHandle hep = he->pair(); // Half-edge pair to he
		
		// if (ver)
		if (verbose || _log_topology_stuff) {
			cout << "n sides before for he and hep:" << he->face()->nsides << " , " << hep->face()->nsides << endl;
		}
		
		VertexHandle v1 = he->from(); // We define v1 and the vertex he points from
		VertexHandle v2 = he->to();	// We define v2 and the vertex he points to

		Vec l = 0.5 * (v2->data().r - v1->data().r); // Vector l points from v1 towards the geometric centre of the v1-v2 line

		Vec rc = v1->data().r + l;

		Vec rot_l = Vec(-l.y, l.x).unit();
		
		Vec new_v1_r = rc - 0.5 * edge_len * rot_l;
		Vec new_v2_r = rc + 0.5 * edge_len * rot_l;
		
		if (verbose || _log_topology_stuff) {
			cout << "Old v1,v2 r : " << v1->data().r << ", " << v2->data().r << endl;
		}
		
		v1->data().r = new_v1_r;
		v2->data().r = new_v2_r;
		
		
		if (verbose || _log_topology_stuff) {
			cout << "New v1,v2 r : " << v1->data().r << ", " << v2->data().r << endl;
		}
		// print()

		v1->set_he(he->idx());
		v2->set_he(hep->idx());

		HEHandle he1 = he->prev();
		HEHandle he2 = he->next();
		HEHandle he3 = hep->prev();
		HEHandle he4 = hep->next();
		
		
		if (verbose || _log_topology_stuff) {
			cout << "Before topo change: " << endl;
			cout << "he1->next(): " << he1->next()->idx() << endl;
			cout << "he2->prev(): " << he2->prev()->idx() << endl;
			cout << "he3->next(): " << he3->next()->idx() << endl;
			cout << "he4->prev(): " << he4->prev()->idx() << endl;
			
			cout << "Target values: " << endl;
			cout << "he2: " << he2->idx() << endl;
			cout << "he1: " << he1->idx() << endl;
			cout << "he4: " << he4->idx() << endl;
			cout << "he3: " << he3->idx() << endl;
		}

		he1->set_next(he2->idx());
		he2->set_prev(he1->idx());
		he3->set_next(he4->idx());
		he4->set_prev(he3->idx());
		
		
		if (verbose || _log_topology_stuff) {
			cout << "After topo change: " << endl;
			cout << "he1->next(): " << he1->next()->idx() << endl;
			cout << "he2->prev(): " << he2->prev()->idx() << endl;
			cout << "he3->next(): " << he3->next()->idx() << endl;
			cout << "he4->prev(): " << he4->prev()->idx() << endl;
		}

		he->set_next(he1->pair()->idx());
		he->set_prev(he4->pair()->idx());
		hep->set_next(he3->pair()->idx());
		hep->set_prev(he2->pair()->idx());
		
		if (verbose || _log_topology_stuff) {
			cout << "Finished setting he, hep next & prev" << endl;
		}
		
		he1->pair()->set_prev(he->idx());
		he2->pair()->set_next(hep->idx());
		he3->pair()->set_prev(hep->idx());
		he4->pair()->set_next(he->idx());
		
		
		if (verbose || _log_topology_stuff) {
			cout << "Finished ssetting he1, he2, he3, he4 pair prev and next" << endl;
		}
		he1->set_to(v2->id);
		he1->pair()->set_from(v2->id);
		
		if (verbose || _log_topology_stuff) {
			cout << "finished setting he1 to & he1 pair from" << endl;
		}
		
		he3->set_to(v1->id);
		he3->pair()->set_from(v1->id);
		
		
		if (verbose || _log_topology_stuff) {
			cout << "finished setting he3 to & he3 pair from" << endl;
		}
		
		he->face()->set_he(he2->idx());
		hep->face()->set_he(he4->idx());
		
		
		if (verbose || _log_topology_stuff) {
			cout << "finihed setting he & hep face he" << endl;
		}
		
		he->set_face(he1->pair()->face()->id);
		hep->set_face(he2->pair()->face()->id);
		
		if (verbose || _log_topology_stuff) {
			cout << "finished setting he & hep face" << endl;
		}
		
		he->face()->nsides = this->face_sides(*(he->face()));
		hep->face()->nsides = this->face_sides(*(hep->face()));
		
		if (verbose || _log_topology_stuff) {
			cout << "finished setting nsides" << endl;
			
			if (he->face()->nsides != 6) {
				cout << "    nsides changed " << endl;
			}
			cout << "  successfully completed a t1 transition!" << endl;
		}

		return true;
	}

	void Mesh::tidyup()
	{
		for (auto& v : _vertices)
			// if (!v.erased)
				v.coordination = this->coordination(v);

		for (auto& e : _edges)
		{
			e.boundary = false;
			if (e.he()->from()->boundary && e.he()->to()->boundary)
				e.boundary = true;
		}
	}

	// Mesh info functions
	double Mesh::area(const Face &f) const
	{
		if (f.outer)
			return 0.0;

		Vec r0 = f.he()->from()->data().r;
		double A = 0.0;
		for (auto he : f.circulator())
		{
			Vec r1 = he.from()->data().r - r0; // this takes care of the boundary conditions
			Vec r2 = he.to()->data().r - r0;
			A += r1.x * r2.y - r2.x * r1.y;
		}
		return 0.5 * fabs(A);
	}

	double Mesh::perim(const Face &f) const
	{
		if (f.outer)
			return 0.0;

		double P = 0.0;
		for (auto he : f.circulator())
        {
			P += (he.to()->data().r - he.from()->data().r).len();
        }

		return P;
	}

	double Mesh::len(const Edge &e)
	{
		return (e.he()->to()->data().r - e.he()->from()->data().r).len();
	}

	int Mesh::coordination(const Vertex &v)
	{
		int i = 0;
		for (auto he : v.circulator()) {
			i++;
		}

		return i;
	}

	int Mesh::face_sides(const Face &f)
	{
		int i = 0;
		for (auto he : f.circulator()) {
			i++;
        }

		return i;
	}

	bool Mesh::is_boundary_face(const Face &f)
	{
		for (auto he : f.circulator())
		{
			if (he.from()->boundary)
				return true;
		}
		return false;
	}

	// // Compute geometric centre of the mesh by tracing positions of boundary vertices
	// Vec Mesh::get_centre()
	// {
	// 	Vec cm(0.0, 0.0);
	// 	for (auto v : _vertices)
	// 		cm += v.data().r;
	// 	return static_cast<double>(1.0 / _vertices.size()) * cm;
	// }

	// // Compute centre of a face
	// Vec Mesh::get_face_centre(const Face &f)
	// {
	// 	if (f.outer)
	// 		return Vec(0.0, 0.0);
	// 	Vec r0 = f.he()->from()->data().r;
	// 	Vec rc(0.0, 0.0);
	// 	for (auto he : f.circulator())
	// 	{
	// 		Vec dr = he.from()->data().r - r0;
	// 		rc += Vec(dr.x, dr.y);
	// 	}
	// 	Vec Rc = (1.0 / f.nsides) * rc + r0;
	// 	return Vec(Rc.x, Rc.y);
	// }

	// Compute position of the face centroid
	Vec Mesh::get_face_centroid(const Face &f)
	{
		if (f.outer)
			return Vec(0.0, 0.0);

		Vec r0 = f.he()->from()->data().r;
		Vec rc(0.0, 0.0);
		for (auto he : f.circulator())
		{
			Vec ri = he.from()->data().r - r0;
			Vec rj = he.to()->data().r - r0;
			double fact = ri.x * rj.y - ri.y * rj.x;
			rc.x += (ri.x + rj.x) * fact;
			rc.y += (ri.y + rj.y) * fact;
		}
		Vec Rc = (1.0 / (6 * this->area(f))) * rc + r0;
		return Vec(Rc.x, Rc.y);
	}

}
#include "types.hpp"
#include "type_circulators.hpp"
#include "mesh.hpp"

namespace VMSim
{
  // Face class ===================================================== 
  HEHandle Face::he() { return _mesh.get_he(_he); }
  HECHandle Face::he() const { return _mesh.get_he(_he); }

  FaceCirculator Face::circulator() { return FaceCirculator(this->he()); }
  FaceCCirculator Face::circulator() const { return FaceCCirculator(this->he()); }

  // HalfEdge class =====================================================
  VertexHandle HalfEdge::from() { return _mesh.get_vertex(_from); }
  VertexHandle HalfEdge::to() { return _mesh.get_vertex(_to); }

  VertexCHandle HalfEdge::from() const { return _mesh.get_vertex(_from); }
  VertexCHandle HalfEdge::to() const { return _mesh.get_vertex(_to); }

  EdgeHandle HalfEdge::edge() { return _mesh.get_edge(_edge); }
  EdgeCHandle HalfEdge::edge() const { return _mesh.get_edge(_edge); }

  FaceHandle HalfEdge::face() { return _mesh.get_face(_face); }
  FaceCHandle HalfEdge::face() const { return _mesh.get_face(_face); }

  HEHandle HalfEdge::pair() { return _mesh.get_he(_pair); }
  HEHandle HalfEdge::next() { return _mesh.get_he(_next); }
  HEHandle HalfEdge::prev() { return _mesh.get_he(_prev); }


  HECHandle HalfEdge::pair() const { return _mesh.get_he(_pair); }
  HECHandle HalfEdge::next() const { return _mesh.get_he(_next); }
  HECHandle HalfEdge::prev() const { return _mesh.get_he(_prev); }
  
  
  Vec HalfEdge::direction() { return this->to()->data().r - this->from()->data().r; }


  // Vertex class ===========================================================
  VertexProperty &Vertex::data() { return _property; }
  VertexProperty Vertex::data() const { return _property; }

  HEHandle Vertex::he() { return _mesh.get_he(_he); }
  HECHandle Vertex::he() const { return _mesh.get_he(_he); }

  VertexCirculator Vertex::circulator() { return VertexCirculator(this->he()); }
  VertexCCirculator Vertex::circulator() const { return VertexCCirculator(this->he()); }

  // Edge class ===========================================================
  HEHandle Edge::he() { return _mesh.get_he(_he); }
  HECHandle Edge::he() const { return _mesh.get_he(_he); }
}
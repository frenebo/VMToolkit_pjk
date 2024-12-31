/*!
 * \file system.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 30-Nov-2023
 * \brief System class 
 */ 

#include "system.hpp"



#include "json.hpp"

using json = nlohmann::json;


using std::cout;
using std::endl;

namespace VMTutorial
{
  void System::read_input_from_jsonstring(const string& json_contents, bool verbose)
  {
    cout << "Reading json string into a json object" << endl;
    if (verbose)
    {
      cout << json_contents << endl;
    }
    
    json j;
    std::istringstream s(json_contents);
    s >> j;
    
    cout << "Finished parsing jon object, now starting to input into the system" << endl;
    
    
    // Check for time step
    if (j["mesh"].find("time_step") != j["mesh"].end())
    {
      if (verbose) cout << "Found time step, setting value from json" << endl;
      _time_step = j["mesh"]["time_step"];
    } else {
      if (verbose) cout << "No time_step found, continuing" << endl;
    }
    
    if (verbose) {
      cout << "Beginning to read vertices" << endl;
    }
    
    // Populate vertices
    for (int i = 0; i < j["mesh"]["vertices"].size(); i++)
    {
      int id = j["mesh"]["vertices"][i]["id"];
      double x = j["mesh"]["vertices"][i]["r"][0];
      double y = j["mesh"]["vertices"][i]["r"][1];
      bool boundary = j["mesh"]["vertices"][i]["boundary"];
      _mesh.add_vertex(Vertex(id, Vec(x, y), boundary, _mesh));
      Vertex &v = _mesh.vertices().back();
      
      // v.erased = j["mesh"]["vertices"][i]["erased"];
      
      if (j["mesh"]["vertices"][i].find("velocity") != j["mesh"]["vertices"][i].end())
      {
        v.data().vel.x = j["mesh"]["vertices"][i]["velocity"][0];
        v.data().vel.y = j["mesh"]["vertices"][i]["velocity"][1];
      }
    }
    cout << "Finished reading vertices." << endl;
    
    // Populate faces
    // bool erased_face = false;
    if (verbose)
    {
      cout << "vertices populated" << endl;
      // log_debug_stats();
      cout << "Beginning to read j -> mesh -> face" << endl;
      cout << "Size of mesh face:" << endl;
      cout << "  " << j["mesh"]["faces"].size() << endl;
    }
    for (int i = 0; i < j["mesh"]["faces"].size(); i++)
    {
      if (verbose) {
        cout << "Reading face at index " << i << endl;
      }
      int face_id = j["mesh"]["faces"][i]["id"];
      
      if (j["mesh"]["faces"][i].find("erased") != j["mesh"]["faces"][i].end())
      {
        // erased_face = j["mesh"]["faces"][i]["erased"];
        // if (erased_face) {
        //   _mesh.add_face(face_id, vector<int>(),  true, verbose);
        // } else {
          const vector<int>& vert_ids = j["mesh"]["faces"][i]["vertices"];
          
          _mesh.add_face(face_id, vert_ids, verbose);
        // }
      }
      else 
      {
        if(verbose) {
          cout << "No erased attribute found, adding face." << endl;
        }
        _mesh.add_face(face_id, j["mesh"]["faces"][i]["vertices"], verbose);
      }
      
      Face& f = _mesh.faces().back();
      
      // if (!erased_face)
      // {
        f.outer = j["mesh"]["faces"][i]["outer"];
      // }
    }
    cout << "Finished reading faces." << endl;
    _mesh.tidyup();
    cout << "Finished mesh setup." << endl;
    cout << "Mesh has " << _mesh.vertices().size() << " vertices " << _mesh.edges().size() << " edges and " << _mesh.faces().size() << " faces." << endl;
    
    _mesh_set = true;
    
    cout << "Finished reading input configuration." << endl;
    
  }
  
  // Python exports
  void export_VertexProperty(py::module& m)
  {
    py::class_<VertexProperty>(m, "VertexProperty")
      .def_readonly("vel", &VertexProperty::vel)
      ;
  }

  void export_EdgeProperty(py::module& m)
  {
    py::class_<EdgeProperty>(m, "EdgeProperty")
      ;
  }

  void export_HEProperty(py::module& m)
  {
    py::class_<HEProperty>(m, "HEProperty")
      ;
  }

  
  void export_FaceProperty(py::module& m)
  {
    py::class_<FaceProperty>(m, "CellProperty")
      ;
  }


  void export_Vertex(py::module& m)
  {
    py::class_<Vertex>(m, "Vertex")
      .def(py::init<Mesh&>())
      .def_readonly("id", &Vertex::id)
      // .def_readonly("erased", &Vertex::erased)
      .def_readwrite("boundary", &Vertex::boundary)
      .def_readonly("coordination", &Vertex::coordination)
      .def("he", [](Vertex& v) { return *(v.he()); })
      .def("vel", [](Vertex& v) { return v.data().vel; })
      .def("property", (VertexProperty& (Vertex::*)()) &Vertex::data, py::return_value_policy::reference);
  }

  void export_VertexCirculator(py::module& m)
  {
    py::class_<VertexCirculator>(m, "VertexCirculator")
        .def(py::init<Vertex>())
        .def("__iter__", &VertexCirculator::__iter__)
        .def("__next__", &VertexCirculator::__next__)
        ;
  }

  void export_Edge(py::module& m)
  {
    py::class_<Edge>(m, "Edge")
      .def(py::init<Mesh&>())
      .def_readonly("i", &Edge::i)
      .def_readonly("j", &Edge::j)
      // .def_readonly("erased", &Edge::erased)
      .def_readonly("boundary", &Edge::boundary)
      .def_property_readonly("id", &Edge::idx)
      .def("he", [](Edge& e) { return *(e.he()); })
      .def("property", (EdgeProperty& (Edge::*)()) &Edge::data, py::return_value_policy::reference)
      ;
  }

  void export_HalfEdge(py::module& m)
  {
    py::class_<HalfEdge>(m, "Junction")
      .def(py::init<Mesh&>())
      .def("vfrom", [](HalfEdge& he) { return *(he.from()); })
      .def("vto", [](HalfEdge& he) { return *(he.to()); })
      .def("edge", [](HalfEdge& he) { return *(he.edge()); })
      .def("id", &HalfEdge::idx)
      .def("next", [](HalfEdge& he) { return *(he.next()); })
      .def("prev", [](HalfEdge& he) { return *(he.prev()); })
      .def("pair", [](HalfEdge& he) { return *(he.pair()); })
      .def("face", [](HalfEdge& he) { return *(he.face()); })
      .def("property", (HEProperty& (HalfEdge::*)()) &HalfEdge::data, py::return_value_policy::reference)
      ;
  }

  void export_Face(py::module& m)
  {
    py::class_<Face>(m, "Cell")
      .def(py::init<Mesh&>())
      .def_readonly("id", &Face::id)
      .def_readonly("neighbours", &Face::nsides)
      .def_readonly("outer", &Face::outer)
      // .def_readonly("erased", &Face::erased)
      .def("he", [](Face& f) { return *(f.he()); })
      .def("property", (FaceProperty& (Face::*)()) &Face::data, py::return_value_policy::reference)
      ;
  }

  void export_FaceCirculator(py::module& m)
  {
    py::class_<FaceCirculator>(m, "FaceCirculator")
        .def(py::init<Face>())
        .def("__iter__", &FaceCirculator::__iter__)
        .def("__next__", &FaceCirculator::__next__)
        ;
  }

  void export_Mesh(py::module& m)
  {
    py::class_<Mesh>(m, "Tissue")
      .def(py::init<>())
      .def("num_vert", &Mesh::num_vert)
      .def("num_cells", &Mesh::num_faces)
      .def("tidyup", &Mesh::tidyup)
      .def("get_vertex", &Mesh::get_vertex , py::return_value_policy::reference)
      .def("get_junction", &Mesh::get_halfedge, py::return_value_policy::reference)
      .def("get_cell", &Mesh::get_face, py::return_value_policy::reference)
      .def("vertices", &Mesh::vertices, py::return_value_policy::reference)
      .def("junctions", &Mesh::edges, py::return_value_policy::reference)
      .def("halfedges", &Mesh::halfedges, py::return_value_policy::reference)
      .def("cells", &Mesh::faces, py::return_value_policy::reference)
      .def("get_cell_centre", [](Mesh& m, int i) -> Vec { return m.get_face_centre(*(m.get_mesh_face(i))); })
      .def("get_cell_centroid", [](Mesh& m, int i) -> Vec { return m.get_face_centroid(*(m.get_mesh_face(i))); })
      .def("get_centre", &Mesh::get_centre)
      .def("cell_area", [](Mesh &m, int i) -> double { return m.area(*(m.get_mesh_face(i))); })
      .def("cell_perim", [](Mesh &m, int i) -> double { return m.perim(*(m.get_mesh_face(i))); })
      .def("get_vertex_positions", &Mesh::get_vertex_positions)
      ;
  }

  void export_System(py::module& m)
  {
    py::class_<System>(m, "System")
      .def(py::init<MyMesh&>())
      .def("read_input_from_jsonstring", &System::read_input_from_jsonstring, py::arg("json_contents"), py::arg("verbose")=false)
      .def("mesh", &System::mesh)
      .def("time_step", &System::time_step)
      .def("simulation_time", &System::simulation_time)
      .def("set_simulation_time_step", &System::set_simulation_time_step)
      .def("log_debug_stats", &System::log_debug_stats)
      ;
  }

  
}


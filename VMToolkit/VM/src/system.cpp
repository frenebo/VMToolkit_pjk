/*!
 * \file system.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 30-Nov-2023
 * \brief System class 
 */ 

#include "system.hpp"

namespace VMTutorial
{
  
  void System::input_from_jsonobj(json& j, bool verbose)
  {
    if (verbose)
    {
      cout << "At start of System::input_from_jsonobj" << endl;
    }
    // // Check if simulation box exists
    // if (j["mesh"].find("box") != j["mesh"].end()) {
    //   if (verbose) cout << "Found box specification" << endl;
      
    // } else {
    //   if (verbose) cout << "no box found in mesh, continuing" << endl;
    // }

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
      _mesh.add_vertex(Vertex<Property>(id, Vec(x, y), boundary, _mesh));
      Vertex<Property> &v = _mesh.vertices().back();
      
      v.erased = j["mesh"]["vertices"][i]["erased"];
      
      if (j["mesh"]["vertices"][i].find("velocity") != j["mesh"]["vertices"][i].end())
      {
        v.data().vel.x = j["mesh"]["vertices"][i]["velocity"][0];
        v.data().vel.y = j["mesh"]["vertices"][i]["velocity"][1];
      }
    }
    cout << "Finished reading vertices." << endl;
    
    // Populate faces
    bool erased_face = false;
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
        erased_face = j["mesh"]["faces"][i]["erased"];
        if (erased_face) {
          _mesh.add_face(face_id, vector<int>(),  true, verbose);
        } else {
          const vector<int>& vert_ids = j["mesh"]["faces"][i]["vertices"];
          
          _mesh.add_face(face_id, vert_ids, false, verbose);
        }
      }
      else 
      {
        if(verbose) {
          cout << "No erased attribute found, adding face." << endl;
        }
        _mesh.add_face(face_id, j["mesh"]["faces"][i]["vertices"], false, verbose);
      }
      
      Face<Property>& f = _mesh.faces().back();
      
      if (!erased_face)
      {
        f.outer = j["mesh"]["faces"][i]["outer"];
      }
    }
    cout << "Finished reading faces." << endl;
    _mesh.tidyup();
    cout << "Finished mesh setup." << endl;
    cout << "Mesh has " << _mesh.vertices().size() << " vertices " << _mesh.edges().size() << " edges and " << _mesh.faces().size() << " faces." << endl;
    
    _mesh_set = true;
    
    cout << "Finished reading input configuration." << endl;
    
  }
  
  void System::read_input_from_jsonstring(const string& json_contents, bool verbose)
  {
    cout << "Reading json string into a json object" << endl;
    if (verbose)
    {
      cout << json_contents << endl;
    }
    json j;
    istringstream s(json_contents);
    s >> j;
    cout << "Finished parsing jon object, now starting to input into the system" << endl;
    input_from_jsonobj(j, verbose);
  }
  
  // Python exports
  void export_VertexProperty(py::module& m)
  {
    py::class_<Property::VertexProperty>(m, "VertexProperty")
      .def_readonly("vel", &Property::VertexProperty::vel)
      ;
  }

  void export_EdgeProperty(py::module& m)
  {
    py::class_<Property::EdgeProperty>(m, "EdgeProperty")
      ;
  }

  void export_HEProperty(py::module& m)
  {
    py::class_<Property::HEProperty>(m, "HEProperty")
      ;
  }

  
  void export_FaceProperty(py::module& m)
  {
    py::class_<Property::FaceProperty>(m, "CellProperty")
      ;
  }


  void export_Vertex(py::module& m)
  {
    py::class_<Vertex<Property>>(m, "Vertex")
      .def(py::init<Mesh<Property>&>())
      .def_readonly("id", &Vertex<Property>::id)
      .def_readonly("erased", &Vertex<Property>::erased)
      .def_readwrite("boundary", &Vertex<Property>::boundary)
      .def_readonly("coordination", &Vertex<Property>::coordination)
      .def("he", [](Vertex<Property>& v) { return *(v.he()); })
      .def("vel", [](Vertex<Property>& v) { return v.data().vel; })
      .def("property", (Property::VertexProperty& (Vertex<Property>::*)()) &Vertex<Property>::data, py::return_value_policy::reference);
  }

  void export_VertexCirculator(py::module& m)
  {
    py::class_<VertexCirculator<Property>>(m, "VertexCirculator")
        .def(py::init<Vertex<Property>>())
        .def("__iter__", &VertexCirculator<Property>::__iter__)
        .def("__next__", &VertexCirculator<Property>::__next__)
        ;
  }

  void export_Edge(py::module& m)
  {
    py::class_<Edge<Property>>(m, "Edge")
      .def(py::init<Mesh<Property>&>())
      .def_readonly("i", &Edge<Property>::i)
      .def_readonly("j", &Edge<Property>::j)
      .def_readonly("erased", &Edge<Property>::erased)
      .def_readonly("boundary", &Edge<Property>::boundary)
      .def_property_readonly("id", &Edge<Property>::idx)
      .def("he", [](Edge<Property>& e) { return *(e.he()); })
      .def("property", (Property::EdgeProperty& (Edge<Property>::*)()) &Edge<Property>::data, py::return_value_policy::reference)
      ;
  }

  void export_HalfEdge(py::module& m)
  {
    py::class_<HalfEdge<Property>>(m, "Junction")
      .def(py::init<Mesh<Property>&>())
      .def("vfrom", [](HalfEdge<Property>& he) { return *(he.from()); })
      .def("vto", [](HalfEdge<Property>& he) { return *(he.to()); })
      .def("edge", [](HalfEdge<Property>& he) { return *(he.edge()); })
      .def("id", &HalfEdge<Property>::idx)
      .def("next", [](HalfEdge<Property>& he) { return *(he.next()); })
      .def("prev", [](HalfEdge<Property>& he) { return *(he.prev()); })
      .def("pair", [](HalfEdge<Property>& he) { return *(he.pair()); })
      .def("face", [](HalfEdge<Property>& he) { return *(he.face()); })
      .def("property", (Property::HEProperty& (HalfEdge<Property>::*)()) &HalfEdge<Property>::data, py::return_value_policy::reference)
      ;
  }

  void export_Face(py::module& m)
  {
    py::class_<Face<Property>>(m, "Cell")
      .def(py::init<Mesh<Property>&>())
      .def_readonly("id", &Face<Property>::id)
      .def_readonly("neighbours", &Face<Property>::nsides)
      .def_readonly("outer", &Face<Property>::outer)
      .def_readonly("erased", &Face<Property>::erased)
      .def("he", [](Face<Property>& f) { return *(f.he()); })
      .def("property", (Property::FaceProperty& (Face<Property>::*)()) &Face<Property>::data, py::return_value_policy::reference)
      ;
  }

  void export_FaceCirculator(py::module& m)
  {
    py::class_<FaceCirculator<Property>>(m, "FaceCirculator")
        .def(py::init<Face<Property>>())
        .def("__iter__", &FaceCirculator<Property>::__iter__)
        .def("__next__", &FaceCirculator<Property>::__next__)
        ;
  }

  void export_Mesh(py::module& m)
  {
    py::class_<Mesh<Property>>(m, "Tissue")
      .def(py::init<>())
      .def("num_vert", &Mesh<Property>::num_vert)
      .def("num_cells", &Mesh<Property>::num_faces)
      .def("tidyup", &Mesh<Property>::tidyup)
      .def("get_vertex", &Mesh<Property>::get_vertex , py::return_value_policy::reference)
      .def("get_junction", &Mesh<Property>::get_halfedge, py::return_value_policy::reference)
      .def("get_cell", &Mesh<Property>::get_face, py::return_value_policy::reference)
      .def("vertices", &Mesh<Property>::vertices, py::return_value_policy::reference)
      .def("junctions", &Mesh<Property>::edges, py::return_value_policy::reference)
      .def("halfedges", &Mesh<Property>::halfedges, py::return_value_policy::reference)
      .def("cells", &Mesh<Property>::faces, py::return_value_policy::reference)
      .def("get_cell_centre", [](Mesh<Property>& m, int i) -> Vec { return m.get_face_centre(*(m.get_mesh_face(i))); })
      .def("get_cell_centroid", [](Mesh<Property>& m, int i) -> Vec { return m.get_face_centroid(*(m.get_mesh_face(i))); })
      .def("get_centre", &Mesh<Property>::get_centre)
      .def("cell_area", [](Mesh<Property> &m, int i) -> double { return m.area(*(m.get_mesh_face(i))); })
      .def("cell_perim", [](Mesh<Property> &m, int i) -> double { return m.perim(*(m.get_mesh_face(i))); })
      .def("get_vertex_positions", &Mesh<Property>::get_vertex_positions)
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


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
    // Check if simulation box exists
    if (j["mesh"].find("box") != j["mesh"].end()) {
      if (verbose) cout << "Found box specification" << endl;
      
      if (j["mesh"]["box"]["periodic"])
      {
        if (j["mesh"]["box"].find("a") != j["mesh"]["box"].end() && j["mesh"]["box"].find("b") != j["mesh"]["box"].end())
        {
          vector<double> a = j["mesh"]["box"]["a"];
          vector<double> b = j["mesh"]["box"]["b"];
          this->set_box(make_shared<Box>(a[0], a[1], b[0], b[1]));
        }
        else if (j["mesh"]["box"].find("lx") != j["mesh"]["box"].end() && j["mesh"]["box"].find("ly") != j["mesh"]["box"].end())
        {
          double lx = j["mesh"]["box"]["lx"];
          double ly = j["mesh"]["box"]["ly"];
          this->set_box(make_shared<Box>(lx, ly));
        }
        else 
          throw runtime_error("There is a problem with the box field in the input JSON file.");
        cout << "Setting periodic simulation box." << endl;
      }
    } else {
      if (verbose) cout << "no box found in mesh, continuing" << endl;
    }

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
      _mesh.add_vertex(Vertex<Property>(id, Vec(x, y, _mesh.box()), boundary, _mesh));
      Vertex<Property> &v = _mesh.vertices().back();
      // this->add_vert_type(j["mesh"]["vertices"][i]["type"]);
      v.erased = j["mesh"]["vertices"][i]["erased"];
      // v.data().vert_type = _vert_types[j["mesh"]["vertices"][i]["type"]];
      // v.data().type_name = get_vert_type_name(v.data().vert_type);
      // v.data().constraint = j["mesh"]["vertices"][i]["constraint"];
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
      log_debug_stats();
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
      // cout << "Reading face index " << i << endl;
      
      // cout 
      if (j["mesh"]["faces"][i].find("erased") != j["mesh"]["faces"][i].end())
      {
        // cout << "ckpt 1" << endl;
        erased_face = j["mesh"]["faces"][i]["erased"];
        if (erased_face) {
          _mesh.add_face(face_id, vector<int>(),  true, verbose);
        } else {
          const vector<int>& vert_ids = j["mesh"]["faces"][i]["vertices"];
		// void add_face(int face_id, const vector<int> &vert_ids, bool erased = false, bool verbose=false);
          _mesh.add_face(face_id, vert_ids, false, verbose);
        }
      }
      else 
      {
        if(verbose) {
          cout << "No erased attribute found, adding face." << endl;
        }
        // cout << "ckpt 2" << endl;
        _mesh.add_face(face_id, j["mesh"]["faces"][i]["vertices"], false, verbose);
      }
      
      // cout << "ckpt 3" << endl;
      Face<Property>& f = _mesh.faces().back();
      f.data().unique_id = f.id;
      
      // cout << "ckpt 4" << endl;
      if (!erased_face)
      {
        // this->add_cell_type(j["mesh"]["faces"][i]["type"]);
        // f.data().face_type = _cell_types[j["mesh"]["faces"][i]["type"]];
        // f.data().type_name = get_cell_type_name(f.data().face_type);
        f.outer = j["mesh"]["faces"][i]["outer"];
        // if (j["mesh"]["faces"][i].find("A0") != j["mesh"]["faces"][i].end())
        //   f.data().A0 = j["mesh"]["faces"][i]["A0"];
        // if (j["mesh"]["faces"][i].find("P0") != j["mesh"]["faces"][i].end())
          // f.data().P0 = j["mesh"]["faces"][i]["P0"];
        if (j["mesh"]["faces"][i].find("n") != j["mesh"]["faces"][i].end())
        {
          double nx = j["mesh"]["faces"][i]["n"][0];
          double ny = j["mesh"]["faces"][i]["n"][1];
          f.data().n = Vec(nx, ny);
        }
        else  
          f.data().n = Vec{0,0};
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
    // if (verboe
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
  
  // Read input from a JSON file
  void System::read_input_from_json_fp(const string& json_file, bool verbose)
  {
    if (verbose) {
      cout << "Reading input json into system" << endl;
      
      log_debug_stats();
    }
    if (_mesh_set)
    {
      cout << "Warning! Mesh has already been set. Overwriting it." << endl;
      _mesh.wipe();
    }
    ifstream inp(json_file.c_str());
    json j;
    inp >> j;
    inp.close();
    
    input_from_jsonobj(j, verbose);
  }


  // Used to be able to make a map of MyoStore
  bool operator<(const VertexHandle<Property>& lv, const VertexHandle<Property>& rv) 
  {
    return (lv->id < rv->id);
  }

  // Python exports

  
  void export_VertexProperty(py::module& m)
  {
    py::class_<Property::VertexProperty>(m, "VertexProperty")
      .def_readonly("vel", &Property::VertexProperty::vel)
      .def_readonly("force", &Property::VertexProperty::force);
  }

  void export_EdgeProperty(py::module& m)
  {
    py::class_<Property::EdgeProperty>(m, "EdgeProperty")
      .def_readwrite("tension", &Property::EdgeProperty::tension);
      // .def_readwrite("l0", &Property::EdgeProperty::l0);
  }

  void export_HEProperty(py::module& m)
  {
    py::class_<Property::HEProperty>(m, "HEProperty")
      .def_readonly("tension", &Property::HEProperty::tension)
      // .def_readonly("l0", &Property::HEProperty::l0)
      .def_readonly("force_type", &Property::HEProperty::force_type);
  }

  
  void export_FaceProperty(py::module& m)
  {
    py::class_<Property::FaceProperty>(m, "CellProperty")
      .def_readonly("unique_id", &Property::FaceProperty::unique_id)
      // .def_readwrite("A0", &Property::FaceProperty::A0)
      // .def_readwrite("P0", &Property::FaceProperty::P0)
      .def_readwrite("n", &Property::FaceProperty::n);
  }


  void export_Vertex(py::module& m)
  {
    py::class_<Vertex<Property>>(m, "Vertex")
      .def(py::init<Mesh<Property>&>())
      .def_readwrite("r", &Vertex<Property>::r)
      .def_readonly("id", &Vertex<Property>::id)
      .def_readonly("erased", &Vertex<Property>::erased)
      .def_readwrite("boundary", &Vertex<Property>::boundary)
      .def_readonly("coordination", &Vertex<Property>::coordination)
      .def("he", [](Vertex<Property>& v) { return *(v.he()); })
      // .def("force", [](Vertex<Property>& v) { return v.data().force; })
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
      .def("property", (Property::EdgeProperty& (Edge<Property>::*)()) &Edge<Property>::data, py::return_value_policy::reference);
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
      .def("property", (Property::HEProperty& (HalfEdge<Property>::*)()) &HalfEdge<Property>::data, py::return_value_policy::reference);
  }

  void export_Face(py::module& m)
  {
    py::class_<Face<Property>>(m, "Cell")
      .def(py::init<Mesh<Property>&>())
      .def_readonly("id", &Face<Property>::id)
      .def_readonly("neighbours", &Face<Property>::nsides)
      .def_readonly("outer", &Face<Property>::outer)
      .def_readonly("erased", &Face<Property>::erased)
      // .def("type", [](Face<Property>& f) { return f.data().face_type; })
      // .def("A0", [](Face<Property>& f) { return f.data().A0; })
      // .def("P0", [](Face<Property>& f) { return f.data().P0; })
      .def("he", [](Face<Property>& f) { return *(f.he()); })
      .def("property", (Property::FaceProperty& (Face<Property>::*)()) &Face<Property>::data, py::return_value_policy::reference);
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
      // .def("set_cell_A0", [](Mesh<Property>& m, int i, double A0) { m.get_face(i).data().A0 = A0;  })
      // .def("set_cell_P0", [](Mesh<Property>& m, int i, double P0) { m.get_face(i).data().P0 = P0;  })
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
      .def("cell_perim", [](Mesh<Property> &m, int i) -> double { return m.perim(*(m.get_mesh_face(i))); });
  }

  void export_System(py::module& m)
  {
    py::class_<System>(m, "System")
      .def(py::init<MyMesh&>())
      .def("read_input_from_json_fp", &System::read_input_from_json_fp, py::arg("input_file"), py::arg("verbose")=false)
      .def("read_input_from_jsonstring", &System::read_input_from_jsonstring, py::arg("json_contents"), py::arg("verbose")=false)
      .def("mesh", &System::mesh)
      .def("time_step", &System::time_step)
      .def("simulation_time", &System::simulation_time)
      .def("set_simulation_time_step", &System::set_simulation_time_step)
      .def("log_debug_stats", &System::log_debug_stats);
  }

  
}


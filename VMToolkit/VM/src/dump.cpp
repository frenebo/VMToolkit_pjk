/*!
 * \file dump.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 30-Nov-2023
 * \brief Dump class
 */

#include "dump.hpp"

namespace VMTutorial
{
	std::string Dump::mesh_to_jsonstr()
	{
		// vector<string> strs = split(mesh_file, '.');
		// string ext = strs[strs.size() - 1];
		// transform(ext.begin(), ext.end(), ext.begin(), [](unsigned char c)
		// 		  { return tolower(c); });
		pt::ptree out;
		pt::ptree rng;
		pt::ptree mesh;
		pt::ptree vertices;
		pt::ptree faces;
		if (_sys.periodic())
		{
			pt::ptree box;
			pt::ptree lx, ly;
			pt::ptree A, B;
			pt::ptree periodic;
			pt::ptree mxx, mxy, myx, myy;
			mxx.put("", _sys.mesh().box()->h._mxx);
			myx.put("", _sys.mesh().box()->h._myx);
			A.push_back(std::make_pair("", mxx));
			A.push_back(std::make_pair("", myx));
			mxy.put("", _sys.mesh().box()->h._mxy);
			myy.put("", _sys.mesh().box()->h._myy);
			B.push_back(std::make_pair("", mxy));
			B.push_back(std::make_pair("", myy));
			lx.put("", _sys.mesh().box()->h._mxx);
			ly.put("", _sys.mesh().box()->h._myy);
			if (fabs(_sys.mesh().box()->h._mxy) >= 1e-6 || fabs(_sys.mesh().box()->h._myx) >= 1e-6)
			{
				box.add_child("a", A);
				box.add_child("b", B);
			}
			else
			{
				box.add_child("lx", lx);
				box.add_child("ly", ly);
			}
			periodic.put("", "true");
			box.add_child("periodic", periodic);
			mesh.add_child("box", box);
		}
		for (auto v : _sys.mesh().vertices())
		{
			pt::ptree vertex;
			pt::ptree id;
			pt::ptree boundary;
			pt::ptree x;
			pt::ptree y;
			// pt::ptree type;
			pt::ptree r;
			pt::ptree fx;
			pt::ptree fy;
			pt::ptree force;
			pt::ptree vx;
			pt::ptree vy;
			pt::ptree velocity;
			pt::ptree constraint;
			pt::ptree erased;
			id.put("", v.id);
			boundary.put("", v.boundary);
			// type.put("", v.data().type_name);
			x.put("", v.r.x);
			y.put("", v.r.y);
			r.push_back(std::make_pair("", x));
			r.push_back(std::make_pair("", y));
			force.push_back(std::make_pair("", fx));
			force.push_back(std::make_pair("", fy));
			vx.put("", v.data().vel.x);
			vy.put("", v.data().vel.y);
			velocity.push_back(std::make_pair("", vx));
			velocity.push_back(std::make_pair("", vy));
			constraint.put("", v.data().constraint);
			erased.put("", v.erased);
			vertex.add_child("id", id);
			vertex.add_child("boundary", boundary);
			// vertex.add_child("type", type);
			vertex.add_child("r", r);
			vertex.add_child("force", force);
			vertex.add_child("velocity", velocity);
			vertex.add_child("constraint", constraint);
			vertex.add_child("erased", erased);
			vertices.push_back(std::make_pair("", vertex));
		}
		mesh.add_child("vertices", vertices);
		int face_id = 0;
		for (auto f : _sys.mesh().faces())
		{
			vector<int> verts;
			vector<double> tension;
			if (!f.erased)
			{
				for (auto he : f.circulator())
				{
					verts.push_back(he.from()->id);
					tension.push_back(he.data().tension);
				} 
			}
			
			pt::ptree face;
			pt::ptree id;
			pt::ptree outer;
			pt::ptree nsides;
			// pt::ptree type;
			pt::ptree erased;
			pt::ptree A0;
			pt::ptree P0;
			pt::ptree fverts;
			pt::ptree rc;
			pt::ptree rc_x, rc_y;

			id.put("", face_id++);
			outer.put("", f.outer);
			nsides.put("", f.nsides);
			// type.put("", f.data().type_name);
			erased.put("", f.erased);
			A0.put("", f.data().A0);
			P0.put("", f.data().P0);
			for (int vid : verts)
			{
				pt::ptree pvid;
				pvid.put("", vid);
				fverts.push_back(std::make_pair("", pvid));
			}
			
			face.add_child("id", id);
			face.add_child("outer", outer);
			// face.add_child("type", type);
			face.add_child("erased", erased);
			face.add_child("A0", A0);
			face.add_child("P0", P0);
			face.add_child("vertices", fverts);
			faces.push_back(std::make_pair("", face));
			//}
		}
		mesh.add_child("faces", faces);

		// save some system properties in the json file
		out.add_child("mesh", mesh);

		// if (ext == "json")
		// {
		std::ostringstream oss;
		pt::write_json(oss, out);
		std::regex reg("\\\"([+-]?(\\d+([.]\\d*)?([eE][+-]?\\d+)?|[.]\\d+([eE][+-]?\\d+)?))\\\"");
		// std::regex reg("\\\"(\\-{0,1}[0-9]+(\\.[0-9]+){0,1})\\\"");
		std::regex reg_bool("\\\"(true|false)\\\"");
		std::string result = std::regex_replace(oss.str(), reg, "$1");
		result = std::regex_replace(result, reg_bool, "$1");
		
		
		return result;
			// std::ofstream file;
			// file.open(mesh_file);
			// file << result;
			// file.close();
		// }
		// else {
			// throw runtime_error("")
		// }
		// else if (ext == "xml")
		// 	pt::write_xml(mesh_file, out, std::locale(), pt::xml_writer_make_settings<string>(' ', 4));
		// else
		// 	pt::write_json(mesh_file, out);
	}


	// // Half-edge
	// void to_json(json &j, const HalfEdge<Property> &he)
	// {
	// 	j = json{{"from", he.from()->id}, {"to", he.to()->id}};
	// }

	// void from_json(const json &j, HalfEdge<Property> &he)
	// {
		
	// }

	// // Edge
	// void to_json(json &j, const Edge<Property> &e)
	// {
	// 	j = json{{"i", e.i}, {"j", e.j}, {"boundary", e.boundary}, {"l0", e.data().l0}};
	// }

	// void from_json(const json &j, Edge<Property> &e)
	// {
	// 	e.i = j.at("i").get<int>();
	// 	e.j = j.at("j").get<int>();
	// 	e.boundary = j.at("boundary").get<bool>();
	// }

	// // Vertex
	// void to_json(json &j, const Vertex<Property> &v)
	// {
	// 	vector<int> neigh;
	// 	for (auto he : v.circulator())
	// 		neigh.push_back(he.to()->id);
		 
	// 	string vert_type;
	// 	j = json{
	// 		{"id", v.id},
	// 		{"r", {v.r.x, v.r.y}},
	// 		{"type", v.data().type_name},
	// 		{"erased", v.erased},
	// 		{"boundary", v.boundary},
	// 		{"constraint", v.data().constraint},
	// 		{"force", {v.data().force.x, v.data().force.y}}
	// 		};
	// }

	// void from_json(const json &j, Vertex<Property> &v)
	// {
	// 	v.id = j.at("id").get<int>();
	// 	v.r.x = j.at("r").get<vector<double>>()[0];
	// 	v.r.y = j.at("r").get<vector<double>>()[1];
	// 	v.data().vert_type = j.at("type").get<int>();
	// 	v.erased = j.at("erased").get<bool>();
	// 	v.boundary = j.at("boundary").get<bool>();
	// 	v.data().constraint = j.at("constraint").get<string>();
	// 	v.data().force.x = j.at("force").get<vector<double>>()[0];
	// 	v.data().force.y = j.at("force").get<vector<double>>()[1];
	// }

	// // Face
	// void to_json(json &j, const Face<Property> &f)
	// {
	// 	vector<int> verts;
	// 	vector<double> lx, ly, ndx, ndy, l0, myo, tension;
	// 	for (const auto he : f.circulator())
	// 	{
	// 		verts.push_back(he.from()->id);
	// 		tension.push_back(he.data().tension);
	// 	}
	// 	j = json{
	// 		{"id", f.id},
	// 		{"outer", f.outer},
	// 		{"nsides", f.nsides},
	// 		{"type", f.data().type_name},
	// 		{"A0", f.data().A0},
	// 		{"P0", f.data().P0},
	// 		{"vertices", verts},
	// 		{"n", {f.data().n.x, f.data().n.y}},
	// 		{"rc", {f.data().rc.x, f.data().rc.y}},
	// 		{"neighbours", f.data().neighs},
	// 		{"tension", tension}};
	// }

	// void from_json(const json &j, Face<Property> &f)
	// {
	// 	f.id = j.at("id").get<int>();
	// 	f.outer = j.at("outer").get<bool>();
	// 	f.nsides = j.at("nsides").get<int>();
	// 	f.data().face_type = j.at("type").get<int>();
	// }

	// // Box
	// void to_json(json &j, const Box &b)
	// {
	// 	vector<double> A = {b.h._mxx, b.h._myx};
	// 	vector<double> B = {b.h._mxy, b.h._myy};
	// 	j = json{
	// 		{"lx", b.h._mxx},
	// 		{"ly", b.h._myy},
	// 		{"a", A},
	// 		{"b", B},
	// 		{"periodic", true}};
	// }

	vector<string> split(const string &s, char delim)
	{
		stringstream ss{s};
		string item;
		vector<string> elems;
		while (std::getline(ss, item, delim))
			elems.push_back(move(item));
		return elems;
	}

	void export_Dump(py::module &m)
	{
		py::class_<Dump>(m, "Dump")
			.def(py::init<System &, ForceCompute &>())
			.def("mesh_to_jsonstr", &Dump::mesh_to_jsonstr)
			// .def("dump_json", &Dump::dump_json)
			.def("set_sfc", &Dump::set_sfc);
	}
}

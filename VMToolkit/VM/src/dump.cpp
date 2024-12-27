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
			// constraint.put("", v.data().constraint);
			erased.put("", v.erased);
			vertex.add_child("id", id);
			vertex.add_child("boundary", boundary);
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
			pt::ptree erased;
			pt::ptree fverts;
			pt::ptree rc;
			pt::ptree rc_x, rc_y;

			id.put("", face_id++);
			outer.put("", f.outer);
			nsides.put("", f.nsides);
			erased.put("", f.erased);
			
			for (int vid : verts)
			{
				pt::ptree pvid;
				pvid.put("", vid);
				fverts.push_back(std::make_pair("", pvid));
			}
			
			face.add_child("id", id);
			face.add_child("outer", outer);
			face.add_child("erased", erased);
			face.add_child("vertices", fverts);
			faces.push_back(std::make_pair("", face));
			//}
		}
		mesh.add_child("faces", faces);

		// save some system properties in the json file
		out.add_child("mesh", mesh);

		std::ostringstream oss;
		pt::write_json(oss, out);
		std::regex reg("\\\"([+-]?(\\d+([.]\\d*)?([eE][+-]?\\d+)?|[.]\\d+([eE][+-]?\\d+)?))\\\"");
		// std::regex reg("\\\"(\\-{0,1}[0-9]+(\\.[0-9]+){0,1})\\\"");
		std::regex reg_bool("\\\"(true|false)\\\"");
		std::string result = std::regex_replace(oss.str(), reg, "$1");
		result = std::regex_replace(result, reg_bool, "$1");
		
		
		return result;
	}
	
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
			.def("set_sfc", &Dump::set_sfc);
	}
}

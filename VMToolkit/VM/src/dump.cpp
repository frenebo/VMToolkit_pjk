/*!
 * \file dump.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 30-Nov-2023
 * \brief Dump class
 */

#include "dump.hpp"
// #include "json.hpp"


#include <string>
#include <sstream>
#include <map>
#include <vector>
#include <regex>

#include <boost/algorithm/string.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/info_parser.hpp>

namespace pt = boost::property_tree;

using std::cerr;
using std::copy;
using std::map;
using std::runtime_error;
using std::vector;

namespace VMTutorial
{
	std::string Dump::mesh_to_jsonstr()
	{
		
		pt::ptree out;
		pt::ptree mesh;
		pt::ptree vertices;
		pt::ptree faces;
		
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
			
			id.put("", v.id);
			boundary.put("", v.boundary);
			x.put("", v.data().r.x);
			y.put("", v.data().r.y);
			r.push_back(std::make_pair("", x));
			r.push_back(std::make_pair("", y));
			force.push_back(std::make_pair("", fx));
			force.push_back(std::make_pair("", fy));
			vx.put("", v.data().vel.x);
			vy.put("", v.data().vel.y);
			velocity.push_back(std::make_pair("", vx));
			velocity.push_back(std::make_pair("", vy));
			
			vertex.add_child("id", id);
			vertex.add_child("boundary", boundary);
			vertex.add_child("r", r);
			vertex.add_child("force", force);
			vertex.add_child("velocity", velocity);
			
			vertices.push_back(std::make_pair("", vertex));
		}
		mesh.add_child("vertices", vertices);
		int face_id = 0;
		for (auto f : _sys.mesh().faces())
		{
			vector<int> verts;
			
			for (auto he : f.circulator())
			{
				verts.push_back(he.from()->id);
			}
			
			pt::ptree face;
			pt::ptree id;
			pt::ptree outer;
			pt::ptree nsides;
			pt::ptree fverts;
			pt::ptree rc;
			pt::ptree rc_x, rc_y;

			id.put("", face_id++);
			outer.put("", f.outer);
			nsides.put("", f.nsides);
			
			for (int vid : verts)
			{
				pt::ptree pvid;
				pvid.put("", vid);
				fverts.push_back(std::make_pair("", pvid));
			}
			
			face.add_child("id", id);
			face.add_child("outer", outer);
			face.add_child("vertices", fverts);
			faces.push_back(std::make_pair("", face));
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
	
	void export_Dump(py::module &m)
	{
		py::class_<Dump>(m, "Dump")
			.def(py::init<System &>())
			.def("mesh_to_jsonstr", &Dump::mesh_to_jsonstr)
			;
	}
}

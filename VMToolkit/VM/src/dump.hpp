/*!
 * \file dump.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 30-Nov-2023
 * \brief Dump class
 */

#ifndef __DUMP_HPP__
#define __DUMP_HPP__

#include <string>
#include <sstream>
#include <fstream>
#include <map>
#include <exception>
#include <iomanip>
#include <memory>
#include <vector>
#include <utility>
#include <algorithm>
#include <cctype>
#include <sys/stat.h>
#include <regex>

#include "json.hpp"

#include "system.hpp"
#include "force_compute.hpp"

using std::back_inserter;
using std::cerr;
using std::copy;
using std::ifstream;
using std::istream_iterator;
using std::istringstream;
using std::map;
using std::move;
using std::ofstream;
using std::runtime_error;
using std::setprecision;
using std::setw;
using std::stod;
using std::stoi;
using std::string;
using std::stringstream;
using std::to_string;
using std::tolower;
using std::transform;
using std::unique_ptr;
using std::vector;

namespace VMTutorial
{

	class Dump
	{
	public:
		Dump(System &sys, ForceCompute &fc) : _sys{sys}, _force_compute{fc}, _sfc(0.95) {}

		// void dump_cells(const string &, bool = false, bool = false);
		// void dump_junctions(const string &, bool = false);
		// void dump_mesh(const string &, bool = false);
		// void dump_json(const string &);
		std::string mesh_to_jsonstr();
		// void dump_cell_directors(const string&, bool = false, bool = false);

		void set_sfc(double sfc) { _sfc = sfc; }

	private:
		System &_sys;
		ForceCompute &_force_compute;
		double _sfc; // Scaling factor for junction output
	};

	vector<string> split(const std::string &, char);

	void export_Dump(py::module &);

}

#endif
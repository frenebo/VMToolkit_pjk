/*!
 * \file dump.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 30-Nov-2023
 * \brief Dump class
 */

#ifndef __DUMP_HPP__
#define __DUMP_HPP__

#include <string>
#include <vector>


#include "system.hpp"


namespace VMTutorial
{

	class Dump
	{
	public:
		Dump(System &sys) : _sys{sys} {}

		std::string mesh_to_jsonstr();

	private:
		System &_sys;
	};

	void export_Dump(py::module &);

}

#endif
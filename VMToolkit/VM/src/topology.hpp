/*!
 * \file topology.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 30-Nov-2023
 * \brief Topology class
 */
 
#ifndef __TOPOLOGY_HPP__
#define __TOPOLOGY_HPP__


#include "system.hpp"
#include "force_compute.hpp"


#include <exception>


using std::exception;
using std::runtime_error;

namespace VMTutorial
{
    class Topology
    {
    public:
        Topology(System &sys, int seed) : _sys{sys},
                                          _min_edge_len{0.02},
                                          _new_edge_len{0.022}
        {
        }
        ~Topology() = default;

        void set_params(const params_type &params)
        {
            for (auto &p : params)
            {
                if (p.first == "min_edge_len") {
                    _min_edge_len = p.second;
                } else if (p.first == "new_edge_len") {
                    _new_edge_len = p.second;
                } else {
                    throw runtime_error("Unknown topoloy flag.");
                }
            }
        };


        void T1();

    private:
        System &_sys;
        double _min_edge_len;
        double _new_edge_len;
    };

    void export_Topology(py::module &);

}


#endif

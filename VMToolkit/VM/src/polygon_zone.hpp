/*!
 * \file polygon_zone.hpp
 * \author Paul Kreymborg, frenebo@gmail.com
 * \date 28-Dec-2024
 * \brief PolygonZone class 
*/

#ifndef __POLYGON_ZONE_HPP__
#define __POLYGON_ZONE_HPP__


#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <string>
#include <vector>

#include "vec.hpp"

using std::vector;
using std::map;

namespace VMTutorial
{
    class PolygonZone
    {
        typedef boost::geometry::model::d2::point_xy<double> BPoint;
        typedef boost::geometry::model::linestring<BPoint> BLineString;
        typedef boost::geometry::model::polygon<BPoint> BPolygon;
        typedef boost::geometry::model::multi_linestring<BLineString> BMultiLinestring;
        
        public:
            PolygonZone() {
                _xmin=0.0;
                _xmax=0.0;
                _ymin=0.0;
                _ymax=0.0;
            }
            
            virtual ~PolygonZone() {}
            
            void set_points(const vector<Vec>& polygon_vertices, bool verbose=false);
            
            double cell_edge_intersection(const Vec& vec_from, const Vec& vec_to,  bool verbose=false) const;
            
            double xmin() const { return _xmin; }
            double xmax() const { return _xmax; }
            double ymin() const { return _ymin; }
            double ymax() const { return _ymax; }
            
        private:
            BPolygon _btpoly;
            double _xmin;
            double _ymin;
            double _xmax;
            double _ymax;
    };
}

#endif
/*!
 * \file polygon_zone.cpp
 * \author Paul Kreymborg, frenebo@gmail.com
 * \date 28-Dec-2024
 * \brief PolygonZone class 
*/

#include "polygon_zone.hpp"

using std::cout;
using std::endl;

namespace VMTutorial
{
    void PolygonZone::set_points(const vector<Vec>& polygon_vertices, bool verbose)
    {
        if (verbose) {
            cout << "PolygonZone::set_points - setting up" << endl;
        }
        
        boost::geometry::clear(_btpoly);
        
        _xmin = polygon_vertices.at(0).x;
        _xmax = _xmin;
        
        _ymin = polygon_vertices.at(0).y;
        _ymax = _ymin;
        
        for (const auto& p_vert_vec : polygon_vertices) {
            double vx = p_vert_vec.x;
            double vy = p_vert_vec.y;
            BPoint vert_pt(vx, vy);
            _btpoly.outer().push_back(vert_pt);
            
            if (verbose) {
                cout << "  add pt " << vx << ", " << vy<< endl;
            }
            
            // if (new_)
            if (vx > _xmax) _xmax = vx;
            if (vx < _xmin) _xmin = vx;
            
            if (vy > _ymax) _ymax = vy;
            if (vy < _ymin) _ymin = vy;
        }
        boost::geometry::correct(_btpoly);
    }
    
    double PolygonZone::cell_edge_intersection(const Vec& vec_from, const Vec& vec_to,  bool verbose) const
    {
        BLineString edge_line;
        edge_line.push_back(BPoint(vec_from.x, vec_from.y));
        edge_line.push_back(BPoint(vec_to.x, vec_to.y));
        
        if (verbose) {
            cout << "          PolygonZone::cell_edge_intersection length of the cell  edge: " << boost::geometry::length(edge_line) << endl;
        }
        
        BMultiLinestring intersection;
        boost::geometry::intersection(_btpoly, edge_line, intersection);
        
        if (verbose) {
            cout << "            finding intersection between _btpoly and edge_line..." << endl;
            cout << "              polygon: ";
            
            //getting the vertices back
            for(auto it = boost::begin(boost::geometry::exterior_ring(_btpoly)); it != boost::end(boost::geometry::exterior_ring(_btpoly)); ++it)
            {
                double x = boost::geometry::get<0>(*it);
                double y = boost::geometry::get<1>(*it);
                cout << " (" << x << "," << y<<") ";
                //use the coordinates...
            }
            cout << endl;
            cout << "              edge_line" << endl;
            
        }
        
        // double edge
        double intersection_length = 0;
        
        for(const auto& intersectionPiece : intersection) {
            // BLineString intersectionPiece = *intersectionIter;
            double piece_intersect_len = boost::geometry::length(intersectionPiece);
            intersection_length += piece_intersect_len;
            if (verbose)  {
                cout << "             found intersection piece, length=" << piece_intersect_len << endl;
            }
            // std::cout << "Piece:" << std::endl;
            // for(auto intersectionPieceIter = intersectionPiece.begin(); intersectionPieceIter != intersectionPiece.end(); ++intersectionPieceIter) {
            // cout << boost::geometry::get<0>(intersectionPiece) << " " << boost::geometry::get<1>(intersectionPiece) << endl;
            // }
        }
        
        return intersection_length;
    }
}

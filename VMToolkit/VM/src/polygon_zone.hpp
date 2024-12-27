#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <string>


#include "vec.hpp"


// namespace bg = boost::geometry;
// namespace bgi = boost::geometry::index;
// namespace

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
            
            void set_points(const vector<Vec>& polygon_vertices, bool verbose=false)
            {
                if (verbose) {
                    cout << "PolygonZone::set_points - setting up" << endl;
                }
                boost::geometry::clear(_btpoly);
                
                // doubl
                // if ()
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
            
            double xmin() const
            {
                return _xmin;
            }
            double xmax() const
            {
                return _xmax;
            }
            double ymin() const
            {
                return _ymin;
            }
            double ymax() const
            {
                return _ymax;
            }
            
            double cell_edge_intersection(const Vec& vec_from, const Vec& vec_to,  bool verbose=false) const
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
            
            virtual ~PolygonZone() {}
        private:
            BPolygon _btpoly;
            double _xmin;
            double _ymin;
            double _xmax;
            double _ymax;
    };
}


#include "force_efield_on_cell_boundary_pixelated.hpp"
#include <set>
#include <cmath>
#include <algorithm>

using std::cout;
using std::endl;
using std::map;
using std::pair;

namespace VMSim
{
    
    void ForceEFieldOnCellBoundPixelated::compute_all_vertex_forces(vector<Vec>& res_forces_out, bool verbose)
    {
        if (! _gridspec) {
            throw runtime_error("ForceEFieldOnCellBoundPixelated::compute_all_vertex_forces - cannot compute forces, the gridspec has not been set yet.");
        }
        
        _clear_compute_cache(verbose);
        _cache_all_computations(verbose);
        
        throw runtime_error("Unimplemented");
    }
    
    Vec ForceEFieldOnCellBoundPixelated::_compute_he_force_for_face(
        const Vertex& vert,
        const HalfEdge& he,
        const Face& f,
        bool verbose=false
    ) {
        throw runtime_error("Unimplemented");
    }
    
    void ForceEFieldOnCellBoundPixelated::_clear_compute_cache(bool verbose)
    {
        if (verbose) {
            cout << "ForceEFieldOnCellBoundPixelated::_clear_compute_cache - clearing the results from llast timestep" << endl;
        }
        
        throw runtime_error("Unimplemented");
    }
    
    void ForceEFieldOnCellBoundPixelated::_cache_all_computations(bool verbose)
    {
        // For each edge
        
        throw runtime_error("Unimplemented");
    }
    
    void ForceEFieldOnCellBoundPixelated::_cache_comp_electric_field_on_edges(bool verbose)
    {
        int n_edges = _sys.cmesh().cedges().size();
        
        set<int> edge_ids_to_calculate;
        
        // Find all the edges belonging to the faces affected by thi force
        for (const auto& face_iter : _face_params) {
            int fid = face_iter.first;
            
            for (const auto& he_it : _sys.cmesh().cfaces().at(fid).circulator()) {
                int edge_idx = he_it.edge().idx();
                
                edge_ids_to_calculate.insert(edge_idx);
            }
        }
        
        
        
        // for (int half)
        
        // for (int edge_idx = 0; edge_idx < n_edges; ++edge_idx) {
        //     const Edge& e = _sys.cmesh().cedges().size().at(edge_idx);
        //     // e.idx();
        // }
        // for (const auto& face : _sys.cmesh().cfaces()) {
        //     _cached_enabled_for_fid.at(fid) = _enabled_for_faceidx(fid, verbose);
        //     _cached_face_perims.at(fid) = _sys.cmesh().perim(face);
            
        //     fid++;
        // }
    }
    
    Vec ForceEFieldOnCellBoundPixelated::_get_pixel_origin_loc(pair<int,int> pix_gridpos)
    {
        double xpos = pix_gridpos.first * _gridspec->grid_spacing() + _grid_spacing->origin_x();
        double ypos = pix_gridpos.second * _gridspec->grid_spacing() + _grid_spacing->origin_y();
        
        return Vec(xpos,ypos);
    }
    
    pair<int, int> ForceEFieldOnCellBoundPixelated::_get_pixel_indices(Vec loc)
    {
        if (! _gridspec) {
            throw runtime_error("ForceEFieldOnCellBoundPixelated::_get_pixel_indices - can't get value, gridspec is unset");
        }
        double xpos_dbl = (loc.x - _gridspec->origin_x()) / _gridspec->grid_spacing();
        double ypos_dbl = loc.y - _gridspec->origin_y() / _gridspec->grid_spacing();
        
        int x_index = static_cast<int>(std::floor(xpos_dbl));
        int y_index = static_cast<int>(std::floor(ypos_dbl));
        
        return pair<int,int>{x_index, y_index};
    }
    
    void ForceEFieldOnCellBoundPixelated::_lineseg_solve_for_y(double x, const Vec& line_start, const Vec& line_end)
    {
        double vx = line_end.x - line_start.x;
        double vy = line_end.y - line_start.y;
        
        double how_far_along_segment = (x - line_start.x) / vx;
        
        return how_far_along_segment * vy + line_start.y;
    }
    
    
    void ForceEFieldOnCellBoundPixelated::_lineseg_solve_for_x(double y, const Vec& line_start, const Vec& line_end)
    {
        double vx = line_end.x - line_start.x;
        double vy = line_end.y - line_start.y;
        
        double how_far_along_segment = (y - line_start.y) / vy;
        
        return how_far_along_segment * vx + line_start.x;
    }
    
    vector<Vec> ForceEFieldOnCellBoundPixelated::_get_pixspec_relative_intersection_positions(const PixIntersectionsSpec& pixspec)
    {
        double pix_size = _gridspec.grid_spacing();
        vector<Vec> positions;
        positions.reserve(4);
        if (pixspec.left_intersection()) {
            positions.push_back(
                Vec(0.0, pixspec.left_intersection())
            );
        }
        if (pixspec.right_intersection()) {
            positions.push_back(
                Vec(pix_size, pixspec.right_intersection())
            );
        }
        if (pixspec.bot_intersection()) {
            positions.push_back(
                Vec(pixsspec.top_intersection(), 0)
            );
        }
        if (pixspec.top_intersection()) {
            positions.push_back(
                Vec(pixspec.top_intersection(), pix_size)
            );
        }
        return positions;
    }
    
    double ForceEFieldOnCellBoundPixelated::_get_intersection_of_line_passthrough_pixel(const PixIntersectionsSpec& pixspec)
    {
        vector<Vec> intersection_rel_positions = _get_pixspec_relative_intersection_positions(pixspec);
        if (intersection_rel_positions.size() != 2) {
            throw runtime_error("Number of intersections must be 2 for a line to pass through a pixel");
        }
        
        Vec intersection_0_rel_pos = intersection_rel_positions.at(0);
        Vec intersection_1_rel_pos = intersection_rel_positions.at(1);
        
        return (intersection_0_rel_pos - intersection_1_rel_pos).len();
    }
    
    double ForceEFieldOnCellBoundPixelated::_get_intersection_of_line_ending_inside_pixel(pair<int,int> pix_gridpos, const PixIntersectionsSpec& pixspec,  const Edge& edge)
    {
        vector<Vec> intersection_rel_positions = _get_pixspec_relative_intersection_positions(pixspec);
        if (intersection_rel_positions.size() != 1) {
            throw runtime_error("Number of intersections must be 1 for a line to go into and end within a pixel");
        }
        
        /*
         * Find which end of this edge is inside the pixel
        */
        Vec edge_start_r = edge.he().from().data().r;
        Vec edge_end_r = edge.he().to().data().r;
        
        pair<int,int> edge_start_gridpos = _get_pixel_indices(edge_start_r);
        pair<int,int> edge_end_gridpos   = _get_pixel_indices(edge_end_r);
        
        Vec contained_endpt_pos(0.0,0.0);
        if (edge_start_gridpos == pix_gridpos) {
            contained_endpt_pos = edge_start_r;
        } else if (edge_end_gridpos == pix_gridpos) {
            contained_endpt_pos = edge_end_r;
        } else {
            throw runtime_error("this is weird - the edge neither starts nor ends in this pixel, but we're in _get_intersection_of_line_ending_inside_pixel ");
        }
        
        /*
         * Calculate the distance between the edge's end point, and the intersection
        */
        Vec intersection_abs_loc = intersection_rel_positions.at(0) + _get_pixel_origin_loc(pix_gridpos);
        
        return (intersection_abs_loc - contained_endpt_pos).len();
    }
    
    Vec ForceEFieldOnCellBoundPixelated::_integrate_field_over_edge(const Edge& edge, bool verbose)
    {
        if (!_field_vals_2d) {
            throw runtime_error("Cannot calculate electric field effect on edge - field value have not been set");
        }
        if (!_gridspec) {
            throw runtime_error("Can't run integrate field over edge - gridspec hasn't been set");
        }
        
        /*
         * First we get all of the pixels this edge passes through, then we calculate the length of edge within
         * each of those pixels. We multiply each sublength by the field at that pixel, to get the result.
         */
        
        map<pair<int,int>, PixIntersectionsSpec> pix_intersections = _get_edge_pixel_intersections(edge, verbose);
        
        // Now that we have all of the pixels that are passed thru...
        int num_cells_with_one_intersections = 0; // This should be exactly two - for the start and end cells
        int num_cells_with_three_intersections = 0; // should be zero
        int num_cells_with_four_intersections = 0; // should be zero
        
        map<pair<int,int>, double> pix_intersected_lengths; // contains the length of edge within each pixel
        
        for (const auto& pit : pix_intersections) {
            const pair<int,int> &       pix_gridpos       = pit->first;
            const PixIntersectionsSpec& pix_intersections = pit->second;
            
            int pix_num_intersections = pix_intersections.count_intersections();
            
            // Keep track of cells with unusual numbers of intersections
            // Zero intersections only happens in one case - when the edge start and ends within one pixel
            if (pix_num_intersections == 0) {
                // this should be the only pixel, if the edge is contained entirely within it
                if (pix_intersections.size() != 1) {
                    cout << "pix_intersections size: " << pix_intersections.size() << endl;
                    throw runtime_error("There is a pixel with no intersections, which isn't the only pixel - error");
                }
                
                double edge_len = (edge.he().to().data().r - edge.he().from().data().r).len();
                
                pix_intersected_lengths[pix_gridpos] = edge_len;
            } else if (pix_num_intersections == 1) {
                num_cells_with_one_intersections++;
                
                pix_intersected_lengths[pix_gridpos] = _get_intersection_of_line_ending_inside_pixel(pix_gridpos, pix_intersections, edge);
            } else if (pix_num_intersections == 2) {
                pix_intersected_lengths[pix_gridpos] = _get_intersection_of_line_passthrough_pixel(pixspec);
            } else {
                cout << " pix n intersections: " << pix_num_intersections << endl;
                throw runtime_error("Unexpected number of intersections...")
            }
        }
        
        Vec tot_summed_field_times_length = Vec(0.0,0.0);
        
        // Sum up field at each pixel, times length of edge within pixel
        for (const auto& it : pix_intersected_lengths) {
            const pair<int,int>& pix_gridpos = it->first;
            double edgelen_within_pixel = it->second;
            
            int pix_pos_x = pix_gridpos->first; // collumn number
            int pix_pos_y = pix_gridpos->second; // row number
            
            // if these pixels lie outside the grid, skip
            if (
                pix_pos_x < 0 ||
                pix_pos_y < 0 ||
                pix_pos_x >= _gridspec->ncells_x() ||
                pix_pos_y >= _gridspec->ncells_y()
            ) {
                continue;
            }
            
            const Vec& field_val = _field_vals_2d->at(pix_pos_x).at(pix_pos_y);
            
            tot_summed_field_times_length += field_val * edgelen_within_pixel;
        }
        
        return tot_summed_field_times_length;
    }
    
    
    
    map<pair<int,int>, PixIntersectionsSpec> ForceEFieldOnCellBoundPixelated::_get_edge_pixel_intersections(
        const Edge& edge, bool verbose,
    ) {
        if (! _gridspec) {
            throw runtime_error("    _get_edge_pixel_intersections - Grid spacing has not been set - can't run");
        }
        
        if (verbose) {
            cout << "ForceEFieldOnCellBoundPixelated::get_edge_pixels - starting" << endl;
        }
        if (! _gridspec) {
            throw runtime_error("    Grid spacing has not been set - can't run");
        }
        
        Vec edge_start_vec = edge.he().from().data().r;
        Vec edge_end_vec = edge.he().to().data().r;
        
        pair<int, int> edge_start_pix_pr = _get_pixel_indices(edge_start_vec);
        int edge_start_x = edge_start_pix_pr.first;
        int edge_start_y = edge_start_pix_pr.second;
        
        pair<int,int> edge_end_pix_pr = _get_pixel_indices(edge_end_vec);
        int edge_end_x = edge_end_pix_pr.first;
        int edge_end_y = edge_end_pix_pr.second;
        
        int edge_min_x = std::min<int>(edge_start_x, edge_end_x);
        int edge_min_y = std::min<int>(edge_start_y, edge_end_y);
        
        int edge_max_x = std::max<int>(edge_start_x, edge_end_x);
        int edge_max_y = std::max<int>(edge_start_y, edge_end_y);
        
        
        // Case 1: the entire edge occupies one pixel
        if (edge_start_x == edge_end_x && edge_start_y == edge_end_y) {
            // This will just be a single pixel, with no intersections
            map<pair<int,int>, PixIntersectionsSpec> just_one_pix_res = {
                {
                    pair<int,int>( edge_start_x, edge_start_y ),
                    PixIntersectionsSpec() // empty intersection spec
                }
            };
            
            return just_one_pix_res;
        }
        // Case 2: the edge will cross at least one boundary between rows or columns:
        //    in THIS case, all pixels that the edge passes through, will be adjacent to these crossings.
        //    Following code calculates crossings and stuff:
        
        /*
         * We start off by finding the columns and rows that will have intersections.
         * If the edge's end points are in separate rows or columns, the edge will have to cross rows or columns somewhere:
         *   For EXAMPLE:
         * If the edge starts in column 1, and ends in column 5, then it must pass through these vertical boundaries
         * between pixel columns: between (1,2), (2,3), (3,4), (4,5)
         * Let's number these intersections by pixels on their right sides...
         * so, the edge that starts in column 1 and ends in column 5, intersects these boundaries:
         * 2,3,4,5
         * Ditto for the rows
         */
        
        // @TODO - what if the edge is nearly vertical, or horizontal... and due to precision error,
        // the cell intersections across column boundaries don't match those for horizontal boundaries?
        // Probablly want to make a note of the problem, then discard cell intersections that don't look right,
        // to avoid crashes. e.g. if a cell has the line passing thru three of its walls it is probably wrong.
        // Or, if the line passees thru one of the walls, but the cell doesn't contain the start or end of the edge - also wrong
        
        vector<int> column_crossing_indices;
        // if edge_min_x = 0, edge_max_x = 3, then column_crossing_indices is 1,2,3
        for (int coli = edge_min_x; coli < edge_max_x; coli++) {
            column_crossing_indices.push_back(coli + 1);
        }
        vector<int> row_crossing_indices;
        for (int rowi = edge_min_y; rowi < edge_max_y; rowi++) {
            row_crossing_indices.push_back(rowi + 1);
        }
        
        if (verbose) {
            cout << "  this edge intersects with these columns: ";
            for (const int& coli : column_crossing_indices) { cout << coli << ","; }
            cout << endl << "    and these rows:";
            for (const int& rowi : row_crossing_indices) { cout << rowi << ","; }
            cout << endl;
        }
        
        /*
            Now that we know the x positions where the edge crosses column boundaries,
            and y positions where the edge crosses row boundaries:
            We need to find the y positions column crossings, and the x positions of row crossings.
        */
        vector<float> column_crossing_y_positions;
        for (const int& col_crossing_idx : column_crossing_indices) {
            // Since crossings are numbered by the column on the rhs of the crossing, the x position
            // of the crosssing is just the x-origin of pixels on the rhs column:
            double col_crossing_x = col_crossing_idx * _gridspec->grid_spacing() + _gridspec->origin_x();
            double col_crossing_y = _lineseg_solve_for_y(col_crossing_x, edge_start_vec, edge_end_vec);
            
            column_crossing_y_positions.push_back(col_crossing_y);
        }
        
        vector<float> row_crossing_x_positions;
        for (const int& row_crossing_idx : row_crossing_indices) {
            
            double row_crossing_y = row_crossing_idx * _gridspec->grid_spacing() + _gridspec->origin_x();
            double row_crossing_x = _lineseg_solve_for_x(row_crossing_y, edge_start_vec, edge_end_vec);
            
            row_crossing_x_positions.push_back(row_crossing_x);
        }
        
        /*
            We now have the positions of all the crossings - so we need to assemble all the cells intersected
        */
        map<pair<int,int>, PixIntersectionsSpec> pix_intersections;
        for (int i=0; i < column_crossing_y_positions.size(); i++) {
            float col_crossing_y_pos = column_crossing_y_positions.at(i);
            int col_crossing_idx = column_crossing_indices.at(i);
            
            int col_crossing_gridrow = static_cast<int>((col_crossing_y_pos - _gridspec->origin_y()) / _gridspec->grid_spacing())
            
            // There are two cells involved in a column crossing:
            // Column crossing index refers to the rhs cell in the crossing, so
            // the columns involved are col_crossing_index - 1, col_crossing_index
            
            pair<int, int> left_cell_gridpos = std::pair<int, int> { col_crossing_gridrow, col_crossing_idx - 1 };
            pair<int, int> right_cell_gridpos = std::pair<int, int> { col_crossing_gridrow, col_crossing_idx };
            
            // Add these cells to the pix intersections if not already present
            if (! pix_intersections.contains(left_cell_gridpos)) {
                pix_intersections[left_cell_gridpos] = PixIntersectionsSpec();
            }
            if (! pix_intersections.contains(right_cell_gridpos)) {
                pix_intersections[right_cell_gridpos] = PixIntersectionsSpec();
            }
            
            // Effectively doing col_crossing_y_position (modulo) grid_spacing
            float crossing_y_position_within_pixel = col_crossing_y_pos - (col_crossing_gridrow * _gridspec->grid_spacing());
            
            // Add the appropriate intersections to each of the cells
            pix_intersections.at(left_cell_gridpos).add_right_intersection(crossing_y_position_within_pixel);
            pix_intersections.at(right_cell_gridpos).add_left_intersection(crossing_y_position_within_pixel);
        }
        
        // Do the same for the row crossings
        for (int i=0; i < row_crossing_x_positions; i++) {
            float row_crossing_x_pos = row_crossing_x_positions.at(i);
            int row_crossing_index = row_crossing_indices.at(i);
            
            int row_crossing_gridcol = static_cast<int>( (row_crossing_x_pos - _gridspec->origin_x()) / _gridspec->grid_spacing() );
            
            pair<int, int> bot_cell_gridpos = std::pair<int,int> { row_crossing_x_pos - 1, row_crossing_gridcol };
            pair<int, int> top_cell_gridpos = std::pair<int,int> { row_crossing_x_pos, row_crossing_gridcol };
            
            if (! pix_intersections.contains(bot_cell_gridpos)) {
                pix_intersections[bot_cell_gridpos] = PixIntersectionsSpec();
            }
            if (! pix_intersections.contains(top_cell_gridpos)) {
                pix_intersections[top_cell_gridpos] = PixIntersectionsSpec();
            }
            
            // Effectively finding row_crossing_x_pos (modulo) grid_spacing
            float crossing_x_position_within_pixel = row_crossing_x_pos - (row_crossing_gridcol * _gridspec->grid_spacing());
            
            pix_intersections.at(bot_cell_gridpos).add_top_intersection(crossing_x_position_within_pixel);
            pix_intersections.at(top_cell_gridpos).add_bot_intersection(crossing_x_position_within_pixel);
        }
        
        return pix_intersections;
    }
    
}
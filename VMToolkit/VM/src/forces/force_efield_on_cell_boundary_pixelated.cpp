

#include "force_efield_on_cell_boundary_pixelated.hpp"
#include <set>
#include <cmath>
#include <algorithm>
#include <set>

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
        
        if (verbose) {
            cout << "  ForceEFieldOnCellBoundPixelated::compute_all_vertex_forces - clearing compute cache" << endl;
        }
        _clear_compute_cache(verbose);
        if (verbose) {
            cout << "  ForceEFieldOnCellBoundPixelated::compute_all_vertex_forces - caching computations" << endl;
        }
        _cache_all_computations(verbose);
        
        res_forces_out.clear();
        int n_vertices = _sys.cmesh().cvertices().size();
        res_forces_out.resize(n_vertices, Vec(0.0, 0.0));
        
        for (const auto& fconf_it : _face_params) {
            int fid = fconf_it.first;
            const EFieldPixFaceConfig& face_config = fconf_it.second;
            const Face& face = _sys.cmesh().cfaces().at(fid);
            
            double face_charge = face_config.charge();
            double face_perimeter = _sys.cmesh().perim(face);
            double charge_density_on_face_perim = face_charge / face_perimeter;
            
            for (const auto& he : face.circulator()) {
                // F = int_{he} dQ_{he} \vec{E(x,y)}
                // With dQ = charge linear density * dl, we get
                // F = (charge linear density) int_{he} dl \vec{E}(x,y)
                // and we already calculated the cached integral, so we just multiply like this:
                
                int edge_idx = he.edge()->idx();
                const Vec& integrated_field_dl = _cached_edge_integrated_field_wrt_length.at(edge_idx);
                
                Vec force_on_half_edge = charge_density_on_face_perim * integrated_field_dl;
                // The force of each half edge is split between each of its vertices
                int vtx_to_id   = he.to()->id;
                int vtx_from_id = he.from()->id;
                
                res_forces_out.at(vtx_to_id)   += force_on_half_edge * 0.5;
                res_forces_out.at(vtx_from_id) += force_on_half_edge * 0.5;
            }
        }
    }
    
    
    void ForceEFieldOnCellBoundPixelated::_clear_compute_cache(bool verbose)
    {
        if (verbose) {
            cout << "ForceEFieldOnCellBoundPixelated::_clear_compute_cache - clearing the results from llast timestep" << endl;
        }
        
        _cached_edge_integrated_field_wrt_length.clear();
    }
    
    void ForceEFieldOnCellBoundPixelated::_cache_all_computations(bool verbose)
    {
        if (verbose) {
            cout << "ForceEFieldOnCellBoundPixelated::_cache_all_computations - starting" << endl;
        }
        if (_cached_edge_integrated_field_wrt_length.size() != 0) {
            throw runtime_error("ForceEFieldOnCellBoundPixelated::_cache_all_computations - cached integrated fieldl wrt length should be empty - was cache cleared before callling this?");
        }
        std::set<int> edges_to_compute;
        
        
        // For each face, find all the edges involved
        for (const auto& fconf_it : _face_params) {
            int fid = fconf_it.first;
            const Face& face = _sys.cmesh().cfaces().at(fid);
            
            for (const auto& he : face.circulator()) {
                int edge_idx = he.edge()->idx();
                
                edges_to_compute.insert(edge_idx);
            }
        }
        
        for (int edge_id : edges_to_compute) {
            const Edge& edge = _sys.cmesh().cedges().at(edge_id);
            Vec integrated_field_for_edge = _integrate_field_over_edge_wrt_length(edge, verbose);
            _cached_edge_integrated_field_wrt_length[edge_id] = integrated_field_for_edge;
        }
    }
    
    
    Vec ForceEFieldOnCellBoundPixelated::_get_pixel_origin_loc(pair<int,int> pix_gridpos) const
    {
        if (!_gridspec) {
            throw runtime_error("ForceEFieldOnCellBoundPixelated::_get_pixel_origin_loc - _gridspec not set");
        }
        double grid_spacing = _gridspec.value().grid_spacing();
        
        double xpos = pix_gridpos.first * grid_spacing + _gridspec.value().origin_x();
        double ypos = pix_gridpos.second * grid_spacing  + _gridspec.value().origin_y();
        
        return Vec(xpos,ypos);
    }
    
    pair<int, int> ForceEFieldOnCellBoundPixelated::_get_pixel_indices(const Vec& loc, bool verbose) const
    {
        if (! _gridspec) {
            throw runtime_error("ForceEFieldOnCellBoundPixelated::_get_pixel_indices - can't get value, gridspec is unset");
        }
        double xpos_dbl = (loc.x - _gridspec.value().origin_x()) / _gridspec.value().grid_spacing();
        double ypos_dbl = ( loc.y - _gridspec.value().origin_y() ) / _gridspec.value().grid_spacing();
        
        int x_index = static_cast<int>(std::floor(xpos_dbl));
        int y_index = static_cast<int>(std::floor(ypos_dbl));
        
        if (verbose) {
            cout << " _get_pixel_indices - x0,y0= " << _gridspec.value().origin_x() << "," << _gridspec.value().origin_y() << endl;
            cout << "   _get_pixel_indices - x,y=" << loc.x <<","  << loc.y << " - so xidx,yidx=" << x_index << "," << y_index << endl;
        }
        
        return pair<int,int>{x_index, y_index};
    }
    
    double ForceEFieldOnCellBoundPixelated::_lineseg_solve_for_y(double x, const Vec& line_start, const Vec& line_end, bool verbose) const
    {
        double vx = line_end.x - line_start.x;
        double vy = line_end.y - line_start.y;
        
        double how_far_along_segment = (x - line_start.x) / vx;
        
        double y_intersect = how_far_along_segment * vy + line_start.y;
        
        
        if (verbose) {
            cout << " _lineseg_solve_for_y: Finding the intersection of the edge with x=" << x << endl;
            cout << "     x0: " << line_start.x << ", y0: " << line_start.y << endl;
            cout << "     vx: " << vx << ", vy: " << vy << endl;
            cout << "     x: " << x << endl;
            cout << "     how far along segment: " << how_far_along_segment << endl;
            cout << "     y intersect: " << y_intersect << endl;
        }
        
        return y_intersect;
    }
    
    
    double ForceEFieldOnCellBoundPixelated::_lineseg_solve_for_x(double y, const Vec& line_start, const Vec& line_end, bool verbose) const
    {
        double vx = line_end.x - line_start.x;
        double vy = line_end.y - line_start.y;
        
        double how_far_along_segment = (y - line_start.y) / vy;
        
        double x_intersect = how_far_along_segment * vx + line_start.x;
        
        if (verbose) {
            cout << " _lineseg_solve_for_x: Finding the intersection of the edge with y=" << y << endl;
            cout << "     x0: " << line_start.x << ", y0: " << line_start.y << endl;
            cout << "     vx: " << vx << ", vy: " << vy << endl;
            cout << "     y: " << y << endl;
            cout << "     how far along segment: " << how_far_along_segment << endl;
            cout << "     x intersect: " << x_intersect << endl;
        }
        
        return x_intersect;
    }
    
    vector<Vec> ForceEFieldOnCellBoundPixelated::_get_pixspec_relative_intersection_positions(const PixIntersectionsSpec& pixspec, bool verbose) const
    {
        if (!_gridspec) {
            throw runtime_error("ForceEFieldOnCellBoundPixelated::_get_pixspec_relative_intersection_positions - cannot calculate when _gridspec has not been set!");
        }
        double pix_size = _gridspec.value().grid_spacing();
        vector<Vec> positions;
        positions.reserve(4);
        if (pixspec.left_intersection()) {
            positions.push_back(
                Vec(0.0, pixspec.left_intersection().value())
            );
        }
        if (pixspec.right_intersection()) {
            positions.push_back(
                Vec(pix_size, pixspec.right_intersection().value())
            );
        }
        if (pixspec.bot_intersection()) {
            positions.push_back(
                Vec(pixspec.bot_intersection().value(), 0)
            );
        }
        if (pixspec.top_intersection()) {
            positions.push_back(
                Vec(pixspec.top_intersection().value(), pix_size)
            );
        }
        return positions;
    }
    
    double ForceEFieldOnCellBoundPixelated::_get_intersection_of_line_passthrough_pixel(const PixIntersectionsSpec& pixspec, bool verbose) const
    {
        if (verbose) {
            cout << "  _get_intersection_of_line_passthrough_pixel - finding intersection" << endl;
        }
        vector<Vec> intersection_rel_positions = _get_pixspec_relative_intersection_positions(pixspec, verbose);
        if (intersection_rel_positions.size() != 2) {
            throw runtime_error("Number of intersections must be 2 for a line to pass through a pixel");
        }
        
        Vec intersection_0_rel_pos = intersection_rel_positions.at(0);
        Vec intersection_1_rel_pos = intersection_rel_positions.at(1);
        
        return (intersection_0_rel_pos - intersection_1_rel_pos).len();
    }
    
    double ForceEFieldOnCellBoundPixelated::_get_intersection_of_line_ending_inside_pixel(
        pair<int,int> pix_gridpos,
        const PixIntersectionsSpec& pixspec, 
        const Edge& edge,
        bool verbose
    ) const {
        if (verbose) {
            cout << "  _get_intersection_of_line_ending_inside_pixel - finding intersection" << endl;
        }
        vector<Vec> intersection_rel_positions = _get_pixspec_relative_intersection_positions(pixspec, verbose);
        if (intersection_rel_positions.size() != 1) {
            throw runtime_error("Number of intersections must be 1 for a line to go into and end within a pixel");
        }
        
        /*
         * Find which end of this edge is inside the pixel
        */
        Vec edge_start_r = edge.he()->from()->data().r;
        Vec edge_end_r = edge.he()->to()->data().r;
        
        pair<int,int> edge_start_gridpos = _get_pixel_indices(edge_start_r, verbose);
        pair<int,int> edge_end_gridpos   = _get_pixel_indices(edge_end_r, verbose);
        
        Vec contained_endpt_pos(0.0,0.0);
        if (edge_start_gridpos == pix_gridpos) {
            contained_endpt_pos = edge_start_r;
        } else if (edge_end_gridpos == pix_gridpos) {
            contained_endpt_pos = edge_end_r;
        } else {
            cout << "_get_intersection_of_line_ending_inside_pixel: " << endl;
            cout << "  Edge start grid pos : " << edge_start_gridpos.first << "," << edge_start_gridpos.second << endl;
            cout << "  Edge end grid pos : " << edge_end_gridpos.first << "," << edge_end_gridpos.second << endl;
            cout << "  Pix grid pos: " << pix_gridpos.first << "," << pix_gridpos.second << endl;
            double relintx = intersection_rel_positions.at(0).x;
            double relinty = intersection_rel_positions.at(0).y;
            cout << "  Intersection relative position within pixel: " << relintx << "," << relinty << endl;
            cout << "  this is weird - the edge neither starts nor ends in this pixel, but we're in _get_intersection_of_line_ending_inside_pixel " << endl;
            cout << "  just returning zero..." << endl;
            return 0.0;
        }
        
        /*
         * Calculate the distance between the edge's end point, and the intersection
        */
        Vec intersection_abs_loc = intersection_rel_positions.at(0) + _get_pixel_origin_loc(pix_gridpos);
        
        return (intersection_abs_loc - contained_endpt_pos).len();
    }
    
    Vec ForceEFieldOnCellBoundPixelated::_integrate_field_over_edge_wrt_length(const Edge& edge, bool verbose) const
    {
        if (verbose) {
            cout << "ForceEFieldOnCellBoundPixelated::_integrate_field_over_edge_wrt_length - starting" << endl;
        }
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
        
        map<pair<int,int>, PixIntersectionsSpec> all_pix_intersections = _get_edge_pixel_intersections(edge, verbose);
        
        // Now that we have all of the pixels that are passed thru...
        int num_cells_with_one_intersections = 0; // This should be exactly two - for the start and end cells
        int num_cells_with_three_intersections = 0; // should be zero
        int num_cells_with_four_intersections = 0; // should be zero
        
        map<pair<int,int>, double> pix_intersected_lengths; // contains the length of edge within each pixel
        
        if (verbose) {
            cout << "Iterating through all pixel intersections... " << endl;
        }
        for (const auto& pit : all_pix_intersections) {
            const pair<int,int> &       pix_gridpos       = pit.first;
            const PixIntersectionsSpec& intersects_for_pix = pit.second;
            
            if (verbose) {
                cout << "  looking at pixel at " << pix_gridpos.first << "," << pix_gridpos.second << endl;
            }
            
            int pix_num_intersections = intersects_for_pix.count_intersections();
            
            // Keep track of cells with unusual numbers of intersections
            // Zero intersections only happens in one case - when the edge start and ends within one pixel
            if (pix_num_intersections == 0) {
                if (verbose) {
                    cout << "  This pixel has no intersections in it!" << endl;
                }
                // this should be the only pixel, if the edge is contained entirely within it
                if (all_pix_intersections.size() != 1) {
                    cout << "pix_intersections size: " << all_pix_intersections.size() << endl;
                    throw runtime_error("There is a pixel with no intersections, which isn't the only pixel - error");
                }
                
                double edge_len = (edge.he()->to()->data().r - edge.he()->from()->data().r).len();
                
                pix_intersected_lengths[pix_gridpos] = edge_len;
            } else if (pix_num_intersections == 1) {
                num_cells_with_one_intersections++;
                
                pix_intersected_lengths[pix_gridpos] = _get_intersection_of_line_ending_inside_pixel(pix_gridpos, intersects_for_pix, edge, verbose);
            } else if (pix_num_intersections == 2) {
                pix_intersected_lengths[pix_gridpos] = _get_intersection_of_line_passthrough_pixel(intersects_for_pix, verbose);
            } else {
                cout << " pix n intersections: " << pix_num_intersections << endl;
                cout << " Unexpected number of intersections! just skipping this one." << endl;
                // throw runtime_error("Unexpected number of intersections...");
            }
        }
        
        Vec tot_summed_field_times_length = Vec(0.0,0.0);
        
        // Sum up field at each pixel, times length of edge within pixel
        for (const auto& it : pix_intersected_lengths) {
            const pair<int,int>& pix_gridpos = it.first;
            double edgelen_within_pixel = it.second;
            
            int pix_pos_x = pix_gridpos.first; // collumn number
            int pix_pos_y = pix_gridpos.second; // row number
            
            // if these pixels lie outside the grid, skip
            if (
                pix_pos_x < 0 ||
                pix_pos_y < 0 ||
                pix_pos_x >= _gridspec.value().ncells_x() ||
                pix_pos_y >= _gridspec.value().ncells_y()
            ) {
                continue;
            }
            
            const Vec& field_val = _field_vals_2d->at(pix_pos_x).at(pix_pos_y);
            
            tot_summed_field_times_length += field_val * edgelen_within_pixel;
        }
        
        return tot_summed_field_times_length;
    }
    
    
    
    map<pair<int,int>, ForceEFieldOnCellBoundPixelated::PixIntersectionsSpec> ForceEFieldOnCellBoundPixelated::_get_edge_pixel_intersections(
        const Edge& edge,
        bool verbose
    ) const {
        if (! _gridspec) {
            throw runtime_error("    _get_edge_pixel_intersections - Grid spacing has not been set - can't run");
        }
        
        if (verbose) {
            cout << "ForceEFieldOnCellBoundPixelated::get_edge_pixels - starting" << endl;
        }
        
        const Vec& edge_start_vec = edge.he()->from()->data().r;
        const Vec& edge_end_vec = edge.he()->to()->data().r;
        
        pair<int, int> edge_start_pix_pr = _get_pixel_indices(edge_start_vec, verbose);
        int edge_start_x = edge_start_pix_pr.first;
        int edge_start_y = edge_start_pix_pr.second;
        
        pair<int,int> edge_end_pix_pr = _get_pixel_indices(edge_end_vec, verbose);
        int edge_end_x = edge_end_pix_pr.first;
        int edge_end_y = edge_end_pix_pr.second;
        
        int edge_min_x = std::min<int>(edge_start_x, edge_end_x);
        int edge_min_y = std::min<int>(edge_start_y, edge_end_y);
        
        int edge_max_x = std::max<int>(edge_start_x, edge_end_x);
        int edge_max_y = std::max<int>(edge_start_y, edge_end_y);
        
        
        // Case 1: the entire edge occupies one pixel
        if (edge_start_x == edge_end_x && edge_start_y == edge_end_y) {
            if (verbose) {
                cout << "ForceEFieldOnCellBoundPixelated::_get_edge_pixel_intersections - it looks like edge starts and ends at the same coordinates - so there is no intersection with any pixel." << endl;
            }
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
            for (const int& coli : column_crossing_indices) {
                cout << coli << " ";
            }
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
            if (verbose) {
                cout << "   Finding the y position of the column crossing at position " << col_crossing_idx << endl;
            }
            // Since crossings are numbered by the column on the rhs of the crossing, the x position
            // of the crosssing is just the x-origin of pixels on the rhs column:
            double col_crossing_x = col_crossing_idx * _gridspec.value().grid_spacing() + _gridspec.value().origin_x();
            double col_crossing_y = _lineseg_solve_for_y(
                col_crossing_x,
                edge_start_vec,
                edge_end_vec,
                verbose
            );
            
            
            column_crossing_y_positions.push_back(col_crossing_y);
        }
        
        vector<float> row_crossing_x_positions;
        for (const int& row_crossing_idx : row_crossing_indices) {
            if (verbose) {
                cout << "  Finding the x position of the row crossing at position " << row_crossing_idx << endl;
            }
            
            double row_crossing_y = row_crossing_idx * _gridspec.value().grid_spacing() + _gridspec.value().origin_y();
            double row_crossing_x = _lineseg_solve_for_x(
                row_crossing_y,
                edge_start_vec,
                edge_end_vec,
                verbose
            );
            
            row_crossing_x_positions.push_back(row_crossing_x);
        }
        
        /*
            We now have the positions of all the crossings - so we need to assemble all the cells intersected
        */
        if (verbose) {
            cout << "Going through column crossings to find pixel intersections" << endl;
            cout << "   origin x=" << _gridspec.value().origin_x() << ", origin y=" << _gridspec.value().origin_y() << ", gridspacing=" << _gridspec.value().grid_spacing() << endl;
        }
        map<pair<int,int>, PixIntersectionsSpec> pix_intersections;
        for (int i=0; i < column_crossing_y_positions.size(); i++) {
            float col_crossing_y_pos = column_crossing_y_positions.at(i);
            int col_crossing_idx = column_crossing_indices.at(i);
            
            double non_floored_col_crossing_gridrow = (col_crossing_y_pos - _gridspec.value().origin_y()) / _gridspec.value().grid_spacing();
            int col_crossing_gridrow = static_cast<int>(std::floor(non_floored_col_crossing_gridrow));
            if (verbose) {
                cout << "     the column crossing at idx #" << col_crossing_idx << " crosses at y=" << col_crossing_y_pos << ", which has gridrow " << col_crossing_gridrow << " (nonfloored=" << non_floored_col_crossing_gridrow << ")" <<endl;
            }
            
            // There are two cells involved in a column crossing:
            // Column crossing index refers to the rhs cell in the crossing, so
            // the columns involved are col_crossing_index - 1, col_crossing_index
            
            pair<int, int> left_cell_gridpos = std::pair<int, int> { col_crossing_idx - 1, col_crossing_gridrow };
            pair<int, int> right_cell_gridpos = std::pair<int, int> { col_crossing_idx, col_crossing_gridrow };
            
            // Add these cells to the pix intersections if not already present
            if (! pix_intersections.contains(left_cell_gridpos)) {
                pix_intersections[left_cell_gridpos] = PixIntersectionsSpec();
            }
            if (! pix_intersections.contains(right_cell_gridpos)) {
                pix_intersections[right_cell_gridpos] = PixIntersectionsSpec();
            }
            
            // Effectively doing col_crossing_y_position (modulo) grid_spacing
            float crossing_y_position_within_pixel = col_crossing_y_pos - (col_crossing_gridrow * _gridspec.value().grid_spacing());
            
            // Add the appropriate intersections to each of the cells
            pix_intersections.at(left_cell_gridpos).add_right_intersection(crossing_y_position_within_pixel);
            pix_intersections.at(right_cell_gridpos).add_left_intersection(crossing_y_position_within_pixel);
        }
        
        // Do the same for the row crossings
        if (verbose) {
            cout << "Going through row crossings to find intersections" << endl;
        }
        for (int i=0; i < row_crossing_x_positions.size(); i++) {
            float row_crossing_x_pos = row_crossing_x_positions.at(i);
            int row_crossing_index = row_crossing_indices.at(i);
            double row_crossing_gridcol_nonfloored = (row_crossing_x_pos - _gridspec.value().origin_x()) / _gridspec.value().grid_spacing();
            int row_crossing_gridcol = static_cast<int>(std::floor( row_crossing_gridcol_nonfloored ));
            if (verbose) {
                cout << "     the row crossing at idx#" << row_crossing_index << " crosses at x=" << row_crossing_x_pos << ", which has gridcol" << row_crossing_gridcol << " (nonfloored=" << row_crossing_gridcol_nonfloored << ")" << endl;
            }
            
            // // @NOTE IMPORTANT!!!
            // // Floating point error makes this occasionally pick the 'wrong cell' - if the intersection is very close to a pixel's corner...
            // double corner_closeness_thresh = 0.0001;
            // // If close to the low end
            // if ((row_crossing_gridcol_nonfloored - row_crossing_gridcol)/_gridspec.value().grid_spacing() < corner_closeness_thresh) {
                
            // }
            // // If close to the high end
            // if ((row_crossing_gridcol_nonfloored - row_crossing_gridcol)/_gridspec.value().grid_spacing() < corner_closeness_thresh) {
                
            // }
            
            
            pair<int, int> bot_cell_gridpos = std::pair<int,int> {
                row_crossing_gridcol,
                row_crossing_index - 1
            };
            pair<int, int> top_cell_gridpos = std::pair<int,int> { 
                row_crossing_gridcol,
                row_crossing_index
            };
            
            if (! pix_intersections.contains(bot_cell_gridpos)) {
                pix_intersections[bot_cell_gridpos] = PixIntersectionsSpec();
            }
            if (! pix_intersections.contains(top_cell_gridpos)) {
                pix_intersections[top_cell_gridpos] = PixIntersectionsSpec();
            }
            
            // Effectively finding row_crossing_x_pos (modulo) grid_spacing
            float crossing_x_position_within_pixel = row_crossing_x_pos - (row_crossing_gridcol * _gridspec.value().grid_spacing());
            
            pix_intersections.at(bot_cell_gridpos).add_top_intersection(crossing_x_position_within_pixel);
            pix_intersections.at(top_cell_gridpos).add_bot_intersection(crossing_x_position_within_pixel);
        }
        
        if (verbose) {
            cout << "Finished finding pixel intersections for this edge!" << endl;
        }
        
        return pix_intersections;
    }
    
}
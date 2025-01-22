#include "new_force_efield_on_cell_boundary_pixelated.hpp"

#include <set>
#include <algorithm>

namespace VMSim
{
namespace PixelatedElectricStuff
{
    using std::set;
    using std::optional;
    
    void NEWForceEFieldOnCellBoundPixelated::compute_all_vertex_forces(vector<Vec>& res, bool verbose)
    {
        _clear_compute_cache(verbose);
        _cache_computations(verbose);
    }
    
    void NEWForceEFieldOnCellBoundPixelated::_clear_compute_cache(bool verbose)
    {
        if (verbose) {
            cout << "NEWForceEFieldOnCellBoundPixelated::_clear_compute_cache - clearing cache" << endl;
        }
        _cached_edge_computation_results.clear();
    }
    
    void NEWForceEFieldOnCellBoundPixelated::_cache_computations(bool verbose)
    {
        if (verbose) {
            cout << "NEWForceEFieldOnCellBoundPixelated::_cache_computations - caching repeatedly needed computations - edge integrations with field" << endl;
        }
        
        set<int> edges_to_compute;
        
        for (const auto& fit : _face_params) {
            int face_id = fit.first;
            
            const Face& face = _sys.cmesh().cfaces().at(face_id);
            
            for (const auto& he : face.circulator()) {
                int edge_idx = he.edge()->idx();
                
                edges_to_compute.insert(edge_idx);
            }
        }
        
        for (int edge_id : edges_to_compute) {
            if (_cached_edge_computation_results.contains(edge_id)) {
                cout << " edge id - " << edge_id << endl;
                throw runtime_error("Error - NEWForceEFieldOnCellBoundPixelated::_cache_computations - trying to cache edge id but it already exists in the cache");
            }
            
            const Edge& edge = _sys.cmesh().cedges().at(edge_id);
            
            _cached_edge_computation_results[edge_id] = _integrate_field_over_edge(edge, verbose);
        }
        
        if (verbose) {
            cout << "NEWForceEFieldOnCellBoundPixelated::_cache_computations - FINISHED" << endl;
        }
    }
    
    vector<double> NEWForceEFieldOnCellBoundPixelated::_get_edge_lengths_passing_thru_pixels(const Edge& edge, const vector<GridCoord>& pixels_intersected_by_edge, bool verbose) const
    {
        if (verbose) {
            cout << "NEWForceEFieldOnCellBoundPixelated::_get_edge_lengths_passing_thru_pixels - finding portions of edge going thru each pixel" << endl;
        }
        /*
         * Finding the segments of the edge that pass through each of these pixels
        */
        vector<double> edge_length_passing_thru_each_pixel;
        
        // Case 1 : the edge starts and ends in the same pixel
        if (pixels_intersected_by_edge.size() == 1) {
            // in this case, the entire edge is contained in the pixel, so we just take its total length:
            Vec edge_start_VEC = edge.he()->from()->data().r;
            Vec edge_end_VEC = edge.he()->to()->data().r;
            double edge_entire_length = (edge_end_VEC - edge_start_VEC).len();
            
            edge_length_passing_thru_each_pixel.push_back(edge_entire_length);
            
            // We're done here, so return
            return edge_length_passing_thru_each_pixel;
        }
        
        // Case 2 : the edge starts and ends in different pixels
        /*
         * Step 1: We want to know where on the pixel the edges pass through them - left, top, right, bottom?
         *    We can tell this by looking at the position of the next pixel in line
        */
        
        // Instantiate array of pixel intersection descriptors, which initialize by default to having no intersections - we'll populate these as we go
        vector<PixelIntersectionsDescriptor> pix_intersections_descriptions(pixels_intersected_by_edge.size(), PixelIntersectionsDescriptor());
        
        for (int pixel_idx = 0; pixel_idx < pixels_intersected_by_edge.size(); pixel_idx++) {
            // if we're looking at the last pixel, it's only intersection will have been with the previous pixel, which already should have been noted
            if (pixel_idx == pixels_intersected_by_edge.size() - 1) {
                if (verbose) {
                    cout << "  Looking at the last pixel in the edge, ending here" << endl;
                }
                continue;
            }
            
            const GridCoord& this_pixel_grid_pos = pixels_intersected_by_edge.at(pixel_idx);
            const GridCoord& next_pixel_grid_pos = pixels_intersected_by_edge.at(pixel_idx + 1);
            
            // If they differ in x
            if (this_pixel_grid_pos.x() != next_pixel_grid_pos.x()  && this_pixel_grid_pos.y() == next_pixel_grid_pos.y()) {
                throw runtime_error("NOT IMPLEMENTED");
            }
            // If they differ in y
            else if (this_pixel_grid_pos.x() == next_pixel_grid_pos.x() && this_pixel_grid_pos.y() != next_pixel_grid_pos.y()) {
                throw runtime_error("NOT IMPLEMENTED");
            }
            // If something weird has happened - normally, consecutive pixels along the edge should only differ in either x or y, not neither or both
            else {
                throw runtime_error("NOT IMPLEMENTED");
            }
        }
        
        if (pixels_intersected_by_edge.size() != edge_length_passing_thru_each_pixel.size()) {
            cout << "pixels_intersected_by_edge.size() - " << pixels_intersected_by_edge.size() << endl;
            cout << "edge_length_passing_thru_each_pixel.size() - " << edge_length_passing_thru_each_pixel.size() << endl;
            throw runtime_error("Lengths don't match - pixels_intersected_by_edge.size() != edge_length_passing_thru_each_pixel.size()");
        }
        
        
        throw runtime_error("Unimplemented");
    }
    
    EdgeComputationResult NEWForceEFieldOnCellBoundPixelated::_integrate_field_over_edge(const Edge& edge, bool verbose) const
    {
        if (verbose) {
            cout << "  NEWForceEFieldOnCellBoundPixelated::_integrate_field_over_edge - starting" << endl;
        }
        vector<GridCoord> pixels_intersected_by_edge = _get_edge_pixel_intersections(edge, verbose);
        
        
        throw runtime_error("Not implemented");
    }
    
    Vec NEWForceEFieldOnCellBoundPixelated::_get_vec_coords_for_grid_pos(GridCoord gc, bool verbose) const
    {
        if (!_gridspec) {
            throw runtime_error("_gridspec not set - _get_vec_coords_for_grid_pos");
        }
        
        int grid_x = gc.x();
        int grid_y = gc.y();
        
        double vec_x = grid_x * _gridspec->spacing_x() + _gridspec->origin_x();
        double vec_y = grid_y * _gridspec->spacing_y() + _gridspec->origin_y();
        
        return Vec(vec_x, vec_y);
    }
    
    GridCoord NEWForceEFieldOnCellBoundPixelated::_get_grid_coords_for_vector(Vec vec, bool verbose) const
    {
        if (!_gridspec) {
            throw runtime_error("_gridspec not set - get_col_row_for_vector");
        }
        
        double rel_x = vec.x - _gridspec->origin_x();
        double rel_y = vec.y - _gridspec->origin_y();
        
        double grid_x_dbl = rel_x / _gridspec->spacing_x();
        double grid_y_dbl = rel_y / _gridspec->spacing_y();
        
        int grid_x = static_cast<int>(std::floor(grid_x_dbl));
        int grid_y = static_cast<int>(std::floor(grid_y_dbl));
        
        return GridCoord(grid_x, grid_y);
    }
    
    vector<int> NEWForceEFieldOnCellBoundPixelated::_get_crossings_generalized(int start_coord, int end_coord, bool verbose) const
    {
        if (start_coord == end_coord) {
            return vector<int>();
        }
        else if (start_coord < end_coord) {
            vector<int> crossings;
            for (int i = start_coord; i < end_coord; i++) {
                crossings.push_back(i + 1);
            }
            return crossings;
        }
        else if (start_coord > end_coord) {
            vector<int> crossings;
            for (int i = start_coord; i > end_coord; i--) {
                crossings.push_back(i);
            }
            return crossings;
        } else {
            throw runtime_error("Something is wrong with these ints");
        }
    }
    
    vector<int> NEWForceEFieldOnCellBoundPixelated::_get_column_crossings_of_edge(GridCoord edge_start, GridCoord edge_end, bool verbose) const
    {
        if (verbose) {
            cout << "  get_column_crossings_of_edge - finding which columns are crossed" << endl;
        }
        
        int start_x = edge_start.x();
        int end_x   = edge_end.x();
        
        return _get_crossings_generalized(start_x, end_x, verbose);
    }
    
    vector<int> NEWForceEFieldOnCellBoundPixelated::_get_row_crossings_of_edge(GridCoord edge_start, GridCoord edge_end, bool verbose) const
    {
        if (verbose) {
            cout << "  get_row_crossings_of_edge - finding which rows are crossed" << endl;
        }
        
        int start_y = edge_start.y();
        int end_y   = edge_end.y();
        
        return _get_crossings_generalized(start_y, end_y, verbose);
    }
    
    double NEWForceEFieldOnCellBoundPixelated::get_relative_position_of_crossing_along_edge(int crossing_coord, double edge_start_coord, double edge_end_coord) const
    {
        return (crossing_coord - edge_start_coord) / (edge_end_coord - edge_start_coord);
    }
    
    vector<GridCoord> NEWForceEFieldOnCellBoundPixelated::_get_edge_pixel_intersections(const Edge& edge, bool verbose) const
    {
        Vec edge_start_VEC = edge.he()->from()->data().r;
        Vec edge_end_VEC = edge.he()->from()->data().r;
        
        GridCoord edge_grid_start = _get_grid_coords_for_vector(edge_start_VEC, verbose);
        GridCoord edge_grid_end = _get_grid_coords_for_vector(edge_end_VEC, verbose);
        
        vector<int> column_crossings_xvals = _get_column_crossings_of_edge(edge_grid_start, edge_grid_end, verbose);
        vector<int> row_crossings_yvals = _get_row_crossings_of_edge(edge_grid_start, edge_grid_end, verbose);
        
        if (verbose) {
            cout << "Edge start: " << edge_grid_start.x() << "," << edge_grid_start.y() << endl;
            cout << "Edge end: " << edge_grid_end.x() << "," << edge_grid_end.y() << endl;
            
            cout << "Column crossings: ";
            for (int i = 0; i < column_crossings_xvals.size(); i++) {
                if (i != 0) cout << ",";
                cout << column_crossings_xvals.at(i);
            }
            cout << endl;
            cout << "Row crossings: ";
            for (int i = 0; i < row_crossings_yvals.size(); i++) {
                if (i != 0) cout << ",";
                cout << row_crossings_yvals.at(i);
            }
            cout << endl;
        }
        // Get relative positions of the row crossings and column crossings, along the edge
        
        // How do we know the starting column of the edge?
        // Case 1: there was a column crossing
        GridCoord current_coord = edge_grid_start;
        vector<double> col_crossing_relative_positions_along_edge;
        for (int col_cross_grid_x : column_crossings_xvals) {
            double col_cross_VEC_x = _get_vec_coords_for_grid_pos(
                GridCoord(col_cross_grid_x, 0),
                verbose
            ).x;
            
            col_crossing_relative_positions_along_edge.push_back(get_relative_position_of_crossing_along_edge(
                col_cross_VEC_x,
                edge_start_VEC.x,
                edge_end_VEC.x
            ));
        }
        vector<double> row_crossing_relative_positions_along_edge;
        for (int row_cross_grid_y : row_crossings_yvals) {
            double row_cross_VEC_y = _get_vec_coords_for_grid_pos(
                GridCoord(0, row_cross_grid_y),
                verbose
            ).y;
            
            row_crossing_relative_positions_along_edge.push_back(get_relative_position_of_crossing_along_edge(
                row_cross_VEC_y,
                edge_start_VEC.y,
                edge_end_VEC.y
            ));
        }
        
        
        
        /* *********************************************
         * Traverse the crossings
         * *********************************************/
        vector<GridCoord> intersected_grid_pixels;
        intersected_grid_pixels.push_back(current_coord);
        
        int col_crossing_index = 0;
        int row_crossing_index = 0;
        while ( col_crossing_index < column_crossings_xvals.size() && row_crossing_index < row_crossings_yvals.size()) {
            // Determine, which crossing is next?
            // case 1 - there are no column crossings left
            //    - next crossing is a row crossing
            // case 2 - there are no row crossings left
            //    - next crossing is a column crossing
            // case 3 - there is at least one col and row crossing left, so we need to calculate which one is next
            
            optional<bool> OPT_next_crossing_is_row;
            if (col_crossing_index == column_crossings_xvals.size()) {
                OPT_next_crossing_is_row = true;
            } else if (row_crossing_index == row_crossings_yvals.size()) {
                OPT_next_crossing_is_row = false;
            } else {
                double col_crossing_pos_on_edge = col_crossing_relative_positions_along_edge.at(col_crossing_index);
                double row_crossing_pos_on_edge = row_crossing_relative_positions_along_edge.at(row_crossing_index);
                
                if (row_crossing_pos_on_edge > col_crossing_pos_on_edge) {
                    OPT_next_crossing_is_row = true;
                } else {
                    OPT_next_crossing_is_row = false;
                }
            }
            
            // check on whether an assignment has happened - shouldn't be necessary but just in case i made am istake
            if (!OPT_next_crossing_is_row) {
                throw runtime_error("hmm???");
            }
            
            bool next_crossing_is_row = OPT_next_crossing_is_row.value();
            
            optional<GridCoord> OPT_next_grid_coord;
            
            if (next_crossing_is_row) {
                int row_crossing_yvalue = row_crossings_yvals.at(row_crossing_index);
                
                // Which direction are we crossing?
                // case 1 - we are currently above the crossing, going downwards
                if (row_crossing_yvalue == current_coord.y()) {
                    OPT_next_grid_coord = GridCoord(
                        current_coord.x(),
                        current_coord.y() - 1
                    );
                }
                // case 2 - we are currently below the crossing, going upwards
                else if (row_crossing_yvalue == current_coord.y() - 1) {
                    OPT_next_grid_coord = GridCoord(
                        current_coord.x(),
                        current_coord.y() + 1
                    );
                } else {
                    cout << "row_crossing_yvalue - " << row_crossing_yvalue << endl;
                    cout << "current_coord.y() - " << current_coord.y() << endl;
                    throw runtime_error("Row crossing doesn't match current stuff");
                }
            } else {
                int col_crossing_xvalue = column_crossings_xvals.at(col_crossing_index);
                
                // Which direction are we crossing
                // case 1 - we are currently to the right of the crossing, going left
                if (col_crossing_xvalue == current_coord.x()) {
                    OPT_next_grid_coord = GridCoord(
                        current_coord.x() - 1,
                        current_coord.y()
                    );
                }
                // case 2 - we are currently to the left of the crossing, going right
                else if (col_crossing_xvalue == current_coord.x() - 1) {
                    OPT_next_grid_coord = GridCoord(
                        current_coord.x() + 1,
                        current_coord.y()
                    );
                } else {
                    cout << "col_crossing_xvalue - " << col_crossing_xvalue << endl;
                    cout << "current_coord.x() - " << current_coord.x() << endl;
                    throw runtime_error("Column crossing doesn't match current stuff");
                }
            }
            
            // check on whether the next grid coord was set - this should always be true, but just checking
            if (!OPT_next_grid_coord) {
                throw runtime_error("OPT_next_grid_coord has not been set ???");
            }
            
            // Increment whichever crossing we decided was next
            if (next_crossing_is_row) {
                row_crossing_index++;
            } else {
                col_crossing_index++;
            }
            
            // Set current grid coord and add to the list
            current_coord = OPT_next_grid_coord.value();
            intersected_grid_pixels.push_back(current_coord);
        }
        
        // We should be at the end
        if (current_coord != edge_grid_end) {
            throw runtime_error("? It seems like we didn't make it to the edge end");
        }
        
        return intersected_grid_pixels;
    }
}
}
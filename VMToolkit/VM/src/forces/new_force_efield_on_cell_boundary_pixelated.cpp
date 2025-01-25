#include "new_force_efield_on_cell_boundary_pixelated.hpp"

#include <set>
#include <algorithm>

namespace VMSim
{
namespace PixelatedElectricStuff
{
    using std::set;
    using std::optional;
    
    void NEWForceEFieldOnCellBoundPixelated::compute_all_vertex_forces(vector<Vec>& vtx_forces_out, bool verbose)
    {
        _clear_compute_cache(verbose);
        _cache_computations(verbose);
        
        int n_vertices = _sys.cmesh().cvertices().size();
        vtx_forces_out.clear();
        vtx_forces_out.resize(n_vertices, Vec(0.0,0.0));
        
        for (const auto& fit : _face_params) {
            int face_id = fit.first;
            const Face& face = _sys.cmesh().cfaces().at(face_id);
            const EFieldPixFaceConfig& faceconf = fit.second;
            
            double face_charge = faceconf.charge();
            double face_perimeter = _sys.cmesh().perim(face);
            
            double perimeter_charge_density = face_charge / face_perimeter;
            
            for (const HalfEdge& he : face.circulator()) {
                int edge_idx = he.edge()->idx();
                const EdgeComputationResult& edge_computation = _cached_edge_computation_results.at(edge_idx);
                
                // force = integral_{over edge} dq * Electric field
                // force = (charge density) * integral_{over edge} dl * electric field
                Vec force_on_half_edge_by_face = perimeter_charge_density * edge_computation.integrated_efield_with_respect_to_dlength();
                
                // Half of the force is applied to each adjacent vertex
                int vertex_idx_from = he.from()->id;
                int vertex_idx_to = he.to()->id;
                
                vtx_forces_out.at(vertex_idx_from) += (1.0/2.0)*force_on_half_edge_by_face;
                vtx_forces_out.at(vertex_idx_to) += (1.0/2.0)*force_on_half_edge_by_face;
            }
        }
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
            
            for (const HalfEdge& he : face.circulator()) {
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
    
    double NEWForceEFieldOnCellBoundPixelated::_get_intersection_with_row_and_find_rel_pos_in_column(const Vec& edge_start_VEC, const Vec& edge_end_VEC, int row_to_intersect, int snap_column, bool verbose) const
    {
        if (!_gridspec) {
            throw runtime_error("_gridspec not set - get_intersection_with_row_and_constrain_to_column can't work");
        }
        GridCoord snap_pixel_gridcoord = GridCoord(snap_column, row_to_intersect);
        Vec pixel_origin = _get_vec_coords_for_grid_pos(snap_pixel_gridcoord, verbose);
        
        double vy = (edge_end_VEC.y - edge_start_VEC.y);
        double vx = (edge_end_VEC.x - edge_start_VEC.x);
        
        double y = pixel_origin.y;
        double dist_along_edge = (y - edge_start_VEC.y) / vy;
        double x_abs = (dist_along_edge * vx) + edge_start_VEC.x;
        
        
        double rel_x = x_abs - pixel_origin.x;
        if (rel_x < 0) {
            rel_x = 0;
        } else if (rel_x > _gridspec->spacing_x()) {
            rel_x = _gridspec->spacing_x();
        }
        
        return rel_x;
    }
    
    double NEWForceEFieldOnCellBoundPixelated::_get_intersection_with_column_and_find_rel_pos_in_row(const Vec& edge_start_VEC, const Vec& edge_end_VEC, int column_to_intersect, int snap_row, bool verbose) const
    {
        if (!_gridspec) {
            throw runtime_error("_gridspec not set - get_intersection_with_column_and_find_rel_pos_in_row can't work");
        }
        GridCoord snap_pixel_gridcoord = GridCoord(column_to_intersect, snap_row);
        Vec pixel_origin = _get_vec_coords_for_grid_pos(snap_pixel_gridcoord, verbose);
        
        double vy = (edge_end_VEC.y - edge_start_VEC.y);
        double vx = (edge_end_VEC.x - edge_start_VEC.x);
        
        double x = pixel_origin.x;
        double dist_along_edge = (x - edge_start_VEC.x) / vx;
        double y_abs = (dist_along_edge * vy) + edge_start_VEC.x;
        
        double rel_y = y_abs - pixel_origin.y;
        
        if (rel_y < 0) {
            rel_y = 0;
        } else if (rel_y > _gridspec->spacing_y()) {
            rel_y = _gridspec->spacing_y();
        }
        
        return rel_y;
    }
    
    vector<Vec> NEWForceEFieldOnCellBoundPixelated::_get_pixel_intersection_absolute_coordinates(const GridCoord& pixel_gridpos, const PixelIntersectionsDescriptor& pix_intersections, bool verbose) const
    {
        vector<Vec> locations_of_intersections;
        Vec pixel_origin = _get_vec_coords_for_grid_pos(pixel_gridpos, verbose);
        
        if (pix_intersections.has_top_intersection()) {
            locations_of_intersections.push_back(Vec(
                pixel_origin.x + pix_intersections.top_intersection_loc(),
                pixel_origin.y + _gridspec->spacing_y()
            ));
        }
        if (pix_intersections.has_bottom_intersection()) {
            locations_of_intersections.push_back(Vec(
                pixel_origin.x + pix_intersections.bottom_intersection_loc(),
                pixel_origin.y
            ));
        }
        if (pix_intersections.has_left_intersection()) {
            locations_of_intersections.push_back(Vec(
                pixel_origin.x,
                pixel_origin.y + pix_intersections.left_intersection_loc()
            ));
        }
        if (pix_intersections.has_right_intersection()) {
            locations_of_intersections.push_back(Vec(
                pixel_origin.x + _gridspec->spacing_x(),
                pixel_origin.y + pix_intersections.right_intersection_loc()
            ));
        }
        
        return locations_of_intersections;
    }
    
    
    vector<double> NEWForceEFieldOnCellBoundPixelated::_get_edge_lengths_passing_thru_pixels(const Edge& edge, const vector<GridCoord>& pixels_intersected_by_edge, bool verbose) const
    {
        if (verbose) {
            cout << "NEWForceEFieldOnCellBoundPixelated::_get_edge_lengths_passing_thru_pixels - finding portions of edge going thru each pixel" << endl;
        }
        /*
         * Finding the segments of the edge that pass through each of these pixels
        */
        
        // Case 1 : the edge starts and ends in the same pixel
        Vec edge_start_VEC = edge.he()->from()->data().r;
        Vec edge_end_VEC = edge.he()->to()->data().r;
        
        if (pixels_intersected_by_edge.size() == 1) {
            // in this case, the entire edge is contained in the pixel, so we just take its total length:
            
            vector<double> result = {
                (edge_end_VEC - edge_start_VEC).len()
            };
            
            // We're done here, so return
            return result;
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
            // get_intersection_with_row_and_snap_to_column
            // If they differ in x
            if (this_pixel_grid_pos.x() != next_pixel_grid_pos.x()  && this_pixel_grid_pos.y() == next_pixel_grid_pos.y()) {
                int shared_row = this_pixel_grid_pos.y();
                
                // Case 1 - the next pixel is to the right
                if (next_pixel_grid_pos.x() == this_pixel_grid_pos.x() + 1) {
                    int intersection_column = next_pixel_grid_pos.x();
                    
                    double intersection_y = _get_intersection_with_column_and_find_rel_pos_in_row(edge_start_VEC, edge_end_VEC, intersection_column, shared_row, verbose);
                    
                    pix_intersections_descriptions.at(pixel_idx).set_right_intersection_loc(intersection_y);
                    pix_intersections_descriptions.at(pixel_idx + 1).set_left_intersection_loc(intersection_y);
                }
                // Case 2 - the next pixel is to the left
                else if (next_pixel_grid_pos.x() == this_pixel_grid_pos.x() - 1) {
                    int intersection_column = this_pixel_grid_pos.x();
                    
                    double intersection_y = _get_intersection_with_column_and_find_rel_pos_in_row(edge_start_VEC, edge_end_VEC, intersection_column, shared_row, verbose);
                    
                    pix_intersections_descriptions.at(pixel_idx).set_left_intersection_loc(intersection_y);
                    pix_intersections_descriptions.at(pixel_idx + 1).set_right_intersection_loc(intersection_y);
                }
                else {
                    throw runtime_error("next pixel is neither directly to the left or right...");
                }
            }
            // If they differ in y
            else if (this_pixel_grid_pos.x() == next_pixel_grid_pos.x() && this_pixel_grid_pos.y() != next_pixel_grid_pos.y()) {
                int shared_column = this_pixel_grid_pos.x();
                
                //Case 1 - the next pixel is above
                if (next_pixel_grid_pos.y() == this_pixel_grid_pos.y() + 1) {
                    int intersection_row = next_pixel_grid_pos.y();
                    
                    double intersection_x = _get_intersection_with_row_and_find_rel_pos_in_column(edge_start_VEC, edge_end_VEC, intersection_row, shared_column, verbose);
                    
                    pix_intersections_descriptions.at(pixel_idx).set_top_intersection_loc(intersection_x);
                    pix_intersections_descriptions.at(pixel_idx + 1).set_bottom_intersection_loc(intersection_x);
                }
                // Case 2 - the next pixel is below
                else if (next_pixel_grid_pos.y() == this_pixel_grid_pos.y() - 1) {
                    int intersection_row = this_pixel_grid_pos.y();
                    
                    double intersection_x = _get_intersection_with_row_and_find_rel_pos_in_column(edge_start_VEC, edge_end_VEC, intersection_row, shared_column, verbose);
                    
                    pix_intersections_descriptions.at(pixel_idx).set_bottom_intersection_loc(intersection_x);
                    pix_intersections_descriptions.at(pixel_idx + 1).set_top_intersection_loc(intersection_x);
                } else {
                    throw runtime_error("next pixel is neither directly above or below...");
                }
            }
            // If something weird has happened - normally, consecutive pixels along the edge should only differ in either x or y, not neither or both
            else {
                cout << "this pixel grid pos: " << this_pixel_grid_pos.x() << "," << this_pixel_grid_pos.y() << endl;
                cout << "next pixel grid pos: " << next_pixel_grid_pos.x() << "," << next_pixel_grid_pos.y() << endl;
                throw runtime_error("Could not determine the neighbor relationship between consecutive pixels alogn edge");
            }
        }
        
        vector<double> edge_length_passing_thru_each_pixel;
        for (int pix_idx = 0; pix_idx < pix_intersections_descriptions.size(); pix_idx++) {
            const PixelIntersectionsDescriptor& pix_intersection_description = pix_intersections_descriptions.at(pix_idx);
            const GridCoord& pix_loc = pixels_intersected_by_edge.at(pix_idx);
            
            // (const GridCoord& pixel_gridpos, const PixelIntersectionsDescriptor& pix_intersections, bool verbose) 
            vector<Vec> pix_intersection_locations = _get_pixel_intersection_absolute_coordinates(pix_loc, pix_intersection_description, verbose);
            // The first pixel
            if (pix_idx == 0) {
                if (pix_intersection_locations.size() != 1) {
                    cout << "pix_intersection_locations.size() - " << pix_intersection_locations.size() << endl;
                    throw runtime_error("Wrong number of pixel intersections for the first pixel - should be exactly one");
                }
                
                Vec intersection_loc = pix_intersection_locations.at(0);
                double len_within_pixel = (edge_start_VEC - intersection_loc).len();
                
                edge_length_passing_thru_each_pixel.push_back(len_within_pixel);
            }
            // The last pixel
            else if (pix_idx == pix_intersections_descriptions.size() - 1) {
                if (pix_intersection_locations.size() != 1) {
                    cout << "pix_intersection_locations.size() - " << pix_intersection_locations.size() << endl;
                    throw runtime_error("Wrong number of pixel intersections for the last pixel - should be exactly one");
                }
                
                Vec intersection_loc = pix_intersection_locations.at(0);
                double len_within_pixel = (edge_end_VEC - intersection_loc).len();
                
                edge_length_passing_thru_each_pixel.push_back(len_within_pixel);
            }
            // Middle pixels
            else {
                if (pix_intersection_locations.size() != 2) {
                    cout << "pix_intersection_locations.size() - " << pix_intersection_locations.size() << endl;
                    throw runtime_error("Wrong number of pixel intersections for middle pixel - should be exactly one");
                }
                
                Vec first_intersection_loc = pix_intersection_locations.at(0);
                Vec second_intersection_loc = pix_intersection_locations.at(1);
                
                double len_within_pixel = (second_intersection_loc - first_intersection_loc).len();
                
                edge_length_passing_thru_each_pixel.push_back(len_within_pixel);
            }
        }
        
        if (pixels_intersected_by_edge.size() != edge_length_passing_thru_each_pixel.size()) {
            cout << "pixels_intersected_by_edge.size() - " << pixels_intersected_by_edge.size() << endl;
            cout << "edge_length_passing_thru_each_pixel.size() - " << edge_length_passing_thru_each_pixel.size() << endl;
            throw runtime_error("Lengths don't match - pixels_intersected_by_edge.size() != edge_length_passing_thru_each_pixel.size()");
        }
        
        return edge_length_passing_thru_each_pixel;
    }
    
    Vec NEWForceEFieldOnCellBoundPixelated::get_electric_field_at_grid_coord(const GridCoord& gc, bool verbose) const
    {
        if (!_gridspec.has_value()) {
            throw runtime_error("get_electric_field_at_grid_coord - can't run, _gridspec has not been set yet");
        }
        if (!_flattened_field_vecs.has_value()) {
            throw runtime_error("get_electric_field_at_grid_coord - can't run, _flattened_field_vecs has not been set yet");
        }
        
        if (gc.x() < 0 || gc.x() >= _gridspec->ncells_x() ||
            gc.y() < 0 || gc.y() >= _gridspec->ncells_y())
        {
            cout << "grid ncells_x:" << _gridspec->ncells_x() << ", ncells_y" << _gridspec->ncells_y() << endl;
            cout << "grid coord: " << gc.x() << "," << gc.y() << endl;
            throw runtime_error("get_electric_field_at_grid_coord - invalid point, outside the grid");
        }
        size_t x = gc.x();
        size_t y = gc.y();
        size_t flattened_pos = x * _gridspec->ncells_y() + y;
        if (flattened_pos >= _flattened_field_vecs->size()) {
            cout << "x:" << x << ", y:" << y << " flattened_pos: " << flattened_pos << endl;
            cout << "_flattened_field_vecs.value().size() - " << _flattened_field_vecs.value().size() << endl;
            throw runtime_error("Flattened pos exceeds size of flattened field vecs");
        }
        
        return _flattened_field_vecs->at(flattened_pos);
    }
    
    EdgeComputationResult NEWForceEFieldOnCellBoundPixelated::_integrate_field_over_edge(const Edge& edge, bool verbose) const
    {
        if (verbose) {
            cout << "  NEWForceEFieldOnCellBoundPixelated::_integrate_field_over_edge - starting" << endl;
        }
        if (!_gridspec) {
            throw runtime_error("NEWForceEFieldOnCellBoundPixelated::_integrate_field_over_edge - err, _gridspec has not been set");
        }
        vector<GridCoord> pixels_intersected_by_edge = _get_edge_pixel_intersections(edge, verbose);
        
        vector<double> edge_lengths_thru_pixels = _get_edge_lengths_passing_thru_pixels(edge, pixels_intersected_by_edge, verbose);
        
        if (edge_lengths_thru_pixels.size() != pixels_intersected_by_edge.size()) {
            throw runtime_error("  Lengths of vectors don't match - edge_lengths_thru_pixels &  pixels_intersected_by_edge");
        }
        
        if (verbose) {
            cout << " pixels intersected by edge: " << endl;
            cout << "gridspec data:" << endl;
            cout << "  ncells_x: " << _gridspec->ncells_x() << endl;
            cout << "  ncells_y: " << _gridspec->ncells_y() << endl;
            cout << "  spacing_x: " << _gridspec->spacing_x() << endl;
            cout << "  spacing_y: " << _gridspec->spacing_y() << endl;
            cout << "  origin_x: " << _gridspec->origin_x() << endl;
            cout << "  origin_y: " << _gridspec->origin_y() << endl;
            
            cout << "  edge start - x,y=" << edge.he()->from()->data().r.x << "," << edge.he()->from()->data().r.y << endl;
            cout << "  edge end - x,y=" << edge.he()->to()->data().r.x << "," << edge.he()->to()->data().r.y << endl;
            
            for (int pidx=0; pidx < pixels_intersected_by_edge.size(); pidx++) {
                const GridCoord& gc = pixels_intersected_by_edge.at(pidx);
                cout << "  pixel x,y,edge length thru pixel - " << gc.x() << "," << gc.y() << "," << edge_lengths_thru_pixels.at(pidx) << endl;
            }
        }
        
        Vec sum_result = Vec(0.0,0.0);
        
        for (int pidx=0; pidx < pixels_intersected_by_edge.size(); pidx++) {
            const GridCoord& gc = pixels_intersected_by_edge.at(pidx);
            // If the pixels lies outside the field, skip it
            if (gc.x() < 0 || gc.x() >= _gridspec->ncells_x() ||
                gc.y() < 0 || gc.y() >= _gridspec->ncells_y()) {
                continue;
            }
            
            double edge_len_within_pix = edge_lengths_thru_pixels.at(pidx);
            
            Vec pix_field = get_electric_field_at_grid_coord(gc, verbose);
            
            sum_result += pix_field * edge_len_within_pix;
        }
        
        return EdgeComputationResult(sum_result);
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
        if (verbose) {
            cout << "_get_crossings_generalized - finding crossing positions" << endl;
        }
        if (start_coord == end_coord) {
            if (verbose) {
                cout << "start and end are the same - returning empty vector" << endl;
            }
            return vector<int>();
        }
        else if (start_coord < end_coord) {
            if (verbose) {
                cout << "start is smaller than end - working forwards" << endl;
            }
            vector<int> crossings;
            for (int i = start_coord; i < end_coord; i++) {
                if (verbose) {
                    cout << "  i:" << i << endl;
                    cout << "     i+1:" << i+1 << endl;
                }
                crossings.push_back(i + 1);
            }
            return crossings;
        }
        else if (start_coord > end_coord) {
            if (verbose) {
                cout << "end is smaller than start - working backwards" << endl;
            }
            vector<int> crossings;
            for (int i = start_coord; i > end_coord; i--) {
                if (verbose) {
                    cout << "   " << i << endl;
                }
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
        if (verbose) {
            cout << "NEWForceEFieldOnCellBoundPixelated::_get_edge_pixel_intersections - starting" << endl;
        }
        Vec edge_start_VEC = edge.he()->from()->data().r;
        Vec edge_end_VEC = edge.he()->to()->data().r;
        
        GridCoord edge_grid_start = _get_grid_coords_for_vector(edge_start_VEC, verbose);
        GridCoord edge_grid_end = _get_grid_coords_for_vector(edge_end_VEC, verbose);
        
        if (edge_grid_start.x() > 1000 || edge_grid_start.y() > 1000) {
            cout << "Edge start vec: " << edge_start_VEC.x << "," << edge_start_VEC.y << endl;
            cout << "Edge end vec: " << edge_end_VEC.x << "," << edge_end_VEC.y << endl;
            throw runtime_error("???");
        }
        
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
        
        // if (verbose) {
        //     cout << "COL CROSSINGS, ROW CROSSINGS ";
        //     for (int col_cross_grid_x : column_crossings_xvals) {
        //         cout << col_cross_grid_x << ",";
        //     }
        //     cout << " - ";
        //     for (int row_cross_grid_y : row_crossings_yvals) {
        //         cout << row_cross_grid_y << ",";
        //     }
        //     cout << endl;
        // }
        
        
        
        /* *********************************************
         * Traverse the crossings
         * *********************************************/
        vector<GridCoord> intersected_grid_pixels;
        intersected_grid_pixels.push_back(current_coord);
        
        int col_crossing_index = 0;
        int row_crossing_index = 0;
        while ( col_crossing_index < column_crossings_xvals.size() || row_crossing_index < row_crossings_yvals.size()) {
            if (verbose) {
                cout << " traversing grid coord - " << current_coord.x() << "," << current_coord.y() << endl;
            }
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
                else if (row_crossing_yvalue == current_coord.y() + 1) {
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
                else if (col_crossing_xvalue == current_coord.x() + 1) {
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
            cout << "CURRENT COORD: " << current_coord.x() << "," << current_coord.y() << endl;
            cout << "EDGE GRID END: " << edge_grid_end.x() << "," << edge_grid_end.y() << endl;
            cout << "EDGE GRID START: " << edge_grid_start.x() << "," << edge_grid_start.y() << endl;
            throw runtime_error("? It seems like we didn't make it to the edge end");
        }
        
        return intersected_grid_pixels;
    }
}
}
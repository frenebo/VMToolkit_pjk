/*!
 * \file force_const_vertex_propulsion.cpp
 * \author Paul Kreymborg, frenebo@gmail.com
 * \date 20-Dec-2024
 * \brief ForceEFieldOnCellBoundary class 
*/ 

#include "force_efield_on_cell_boundary_uniform.hpp"

#include <limits> // std::numeric_limits<double>::signaling_NaN

using std::runtime_error;
using std::cout;
using std::endl;


namespace VMSim
{
    
	void ForceEFieldOnCellBoundary::compute_all_vertex_forces(vector<Vec>& res, bool verbose)
	{
		if (verbose) {
			cout << "   ForceEFieldOnCellBoundary::compute_all_vertex_forces - starting" << endl;
        }
        
		_clear_compute_cache(verbose);
		_cache_mesh_computations(verbose);
		
		size_t n_vertices = _sys.cmesh().cvertices().size();	
		res.resize(n_vertices, Vec(0.0,0.0));
			
		for (const auto& vertex : _sys.cmesh().cvertices())
		{
            Vec v_force(0.0,0.0);
			for (const auto& he : vertex.circulator()) {
				if (verbose) {
					cout << "    computing force on vertex " << vertex.id << " by halfedge " << he.idx() << endl;
				}
				
				v_force += _compute_he_force(vertex, he, verbose);
			}
            res.at(vertex.id) = v_force;
		}
	}
    
    void ForceEFieldOnCellBoundary::_clear_compute_cache(bool verbose)
    {
        if (verbose) { cout << "     ForceEFieldOnCellBoundary::_clear_compute_cache - clearing cache" << endl; }
        _cached_vertices_outside_bbox.clear();
        _cached_enabled_for_fid.clear();
        _cached_face_perims.clear();
        _cached_edgelength_intersecting_polygon_by_vtx_ids.clear();
    }
    
    void ForceEFieldOnCellBoundary::_cache_mesh_computations(bool verbose)
    {        
        if (verbose) { cout << "     ForceEFieldOnCellBoundary::_cache_mesh_computations - precalculating some reused values to speed stuff up " << endl; }

        if (verbose) { cout << "         caching vertex values" << endl; }
        size_t n_vertices = _sys.cmesh().cvertices().size();
        _cached_vertices_outside_bbox.resize(n_vertices, false);
        size_t vid = 0;
        for (const auto& vtx : _sys.cmesh().cvertices()) {
            
            // cout
            _cached_vertices_outside_bbox.at(vid) = _vec_outside_polygon_bbox(vtx.data().r, verbose);
            
            vid++;
        }
        
        if (verbose) { cout << "         caching face values" << endl; }
        // Calculate stuff relating to the faces
        size_t n_faces = _sys.cmesh().cfaces().size();
        
        _cached_enabled_for_fid.resize(n_faces, false);
        
        double sNaN = std::numeric_limits<double>::signaling_NaN();
        _cached_face_perims.resize(n_faces, sNaN);
        
        size_t fid = 0;
        for (const auto& face : _sys.cmesh().cfaces()) {
            _cached_enabled_for_fid.at(fid) = _enabled_for_faceidx(fid, verbose);
            _cached_face_perims.at(fid) = _sys.cmesh().perim(face);
            
            fid++;
        }
    }
    
    Vec ForceEFieldOnCellBoundary::_compute_he_force_for_face(
        const Vertex& vert,
        const HalfEdge& he,
        const Face& f,
        double edge_len_within_polygon,
        bool verbose
    ) {
        if (verbose) {
            cout << "          force enabled for fid " << f.id << ", calculating electric force" << endl;
        }
        double Q1 = _cell_charges_by_face_index.at(f.id);
        double P1 = _cached_face_perims.at(f.id);
        double F1_charge_within_poly = Q1 * (edge_len_within_polygon / P1);
        
        if (verbose) {
            cout << "          edge length inside polygon: " << edge_len_within_polygon << " total perim: " << P1 << " charge within poly " << Q1 << endl;
        }
        
        return Vec(
            F1_charge_within_poly * _E_x_param,
            F1_charge_within_poly * _E_y_param
        );
    }
    
    void ForceEFieldOnCellBoundary::set_global_params(const params_type& num_params, const map<string,string>& str_params, bool verbose)
    {
        string field_type = str_params.at("field_type");
        if (field_type != "constant") {
            throw runtime_error("Unknown field type '" + field_type + "'");
        }
        
        double E_x = num_params.at("E_x");
        double E_y = num_params.at("E_y");
        
        string region_type = str_params.at("region_type");
        if (region_type != "polygon") {
            throw runtime_error("Unknown region type '" + region_type + "'");
        }
        
        vector<Vec> polygon_vertices;
        
        int n_polygon_vertices = num_params.at("n_polygon_vertices");
        for (int vert_i = 0; vert_i < n_polygon_vertices; ++vert_i) {
            double pvert_x = num_params.at("poly_x" + std::to_string(vert_i));
            double pvert_y = num_params.at("poly_y" + std::to_string(vert_i));
            polygon_vertices.push_back(Vec(pvert_x, pvert_y));
        }
        
        _E_y_param = E_y;
        _E_x_param = E_x;
        
        _poly_zone.set_points(polygon_vertices, verbose);
    }
    
    
    
    bool ForceEFieldOnCellBoundary::_vec_outside_polygon_bbox(const Vec& v, bool verbose)
    {
        bool vec_inside_bounding_box = (
            v.x <= _poly_zone.xmax() &&
            v.x >= _poly_zone.xmin() &&
            v.y <= _poly_zone.ymax() &&
            v.y >= _poly_zone.ymin()
        );
        
        return ! vec_inside_bounding_box;
    }
    
    bool ForceEFieldOnCellBoundary::_enabled_for_faceidx(int fid, bool verbose)
    {
        if (fid >= _force_enabled_mask_by_face_index.size()) return false;
        
        bool is_enabled = _force_enabled_mask_by_face_index[fid];
        
        if (verbose && is_enabled) {
            cout << "            Confirmed that efield force is on for fid=" << fid << endl;
        }
        
        return is_enabled;
    }
    
    void ForceEFieldOnCellBoundary::set_face_params_facewise(const vector<int>& fids, const vector<params_type>& params, bool verbose)
    {
        int max_fid = 0;
        for (const auto& fid : fids) {
            if (fid > max_fid) max_fid = fid;
        }
        
        if (max_fid >= _force_enabled_mask_by_face_index.size()) {
            _force_enabled_mask_by_face_index.resize(max_fid + 1, false);
            _cell_charges_by_face_index.resize(max_fid + 1, 0.0);
        }
        
        for (size_t i = 0; i < fids.size(); ++i) {
            int fid = fids.at(i);
            const params_type& fparam = params.at(i);
            
            _force_enabled_mask_by_face_index.at(fid) = true;
            _cell_charges_by_face_index.at(fid) = 1.0;
        }
    }
    
    double ForceEFieldOnCellBoundary::_lazy_load_cell_edge_intersection(int vid_from, int vid_to, bool verbose) {
        std::pair<int, int> edge_vid_pair(vid_from, vid_to);
        
        map<std::pair<int,int>,double>::iterator it_cached = _cached_edgelength_intersecting_polygon_by_vtx_ids.find(edge_vid_pair);
        // If this edge has been cached, return
        if(it_cached != _cached_edgelength_intersecting_polygon_by_vtx_ids.end())
        {
            return it_cached->second;
        }
        
        const Vec& v_r = _sys.cmesh().cvertices().at(vid_from).data().r;
        const Vec& vto_r = _sys.cmesh().cvertices().at(vid_to).data().r;
        
        // otherwise calculate the intersection to remember, then return value
        double edge_len_within_polygon = _poly_zone.cell_edge_intersection(v_r, vto_r, verbose);
        
        std::pair<int, int> rev_edge_vid_pair(vid_to, vid_from);
        _cached_edgelength_intersecting_polygon_by_vtx_ids[edge_vid_pair] = edge_len_within_polygon;
        _cached_edgelength_intersecting_polygon_by_vtx_ids[rev_edge_vid_pair] = edge_len_within_polygon;
        
        return edge_len_within_polygon;
    }
     
    Vec ForceEFieldOnCellBoundary::_compute_he_force(const Vertex& vert, const HalfEdge& he, bool verbose)
    {
        if (verbose) {
            cout << "        ForceEFieldOnCellBoundary::_compute_he_force - finding force by he " << he.idx() << endl;
        }
        
        int vfrom_id = vert.id;
        int vto_id = he.to()->id;
        
        Vec v_r = vert.data().r;
        Vec vto_r = he.to()->data().r;
        
        // If both ends of this half edge lie outside the rectangle containing the electric field,
        // then the edge cannot possibly be experiencing any force.
        if (
            _cached_vertices_outside_bbox[vfrom_id] &&
            _cached_vertices_outside_bbox[vto_id]
        ) {
            return Vec(0.0, 0.0);
        }
        
        double edge_len_within_polygon = _lazy_load_cell_edge_intersection(vfrom_id, vto_id, verbose);
        
        const Face& f   = *(he.face());         // cell to the right of the half edge
        const Face& fp  = *(he.pair()->face()); // pair cell (opposite side of the same junction)
        
        // Total force on vertex by the half edge, x and y
        Vec hev_force(0.0,0.0);
        
        if (_cached_enabled_for_fid.at(f.id)) {
            if (verbose) {
                cout << "          force enabled for fid " << f.id << ", calculating electric force" << endl;
            }
            
            hev_force += _compute_he_force_for_face(vert, he, f, edge_len_within_polygon, verbose);
        } else {
            if (verbose) { cout << "          force not enabled for fid " << f.id << ", skipping..." << endl; }
        }
        if (_cached_enabled_for_fid.at(fp.id)) {
            if (verbose) { cout << "          force enabled for fid " << fp.id << ", calculating electric force" << endl; }
            
            hev_force += _compute_he_force_for_face(vert, he, fp, edge_len_within_polygon, verbose);
        } else {
            if (verbose) { cout << "          force not enabled for fid " << fp.id << ", skipping..." << endl; }
        }
        
        if (verbose) {
            cout << "        ForceEFieldOnCellBoundary::_compute_he_force - force on vertex " << vert.id << " by he " << he.idx() << " is :" << hev_force.x << ", " << hev_force.y << endl;
        }
        
        return hev_force;
    }
}
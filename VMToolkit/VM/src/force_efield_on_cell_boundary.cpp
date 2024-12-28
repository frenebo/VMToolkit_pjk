
#include "force_efield_on_cell_boundary.hpp"


namespace VMTutorial
{
    Vec ForceEFieldOnCellBoundary::compute_he_force_for_face(
        const Vertex<Property>& vert,
        const HalfEdge<Property>& he,
        const Face<Property>& f,
        double edge_len_within_polygon,
        bool verbose
    ) {
        if (verbose) {
            cout << "          force enabled for fid " << f.id << ", calculating electric force" << endl;
        }
        double Q1 = _cell_charges_by_face_index.at(f.id);
        double P1 = _sys.mesh().perim(f);
        double F1_charge_within_poly = Q1 * (edge_len_within_polygon / P1);
        if (verbose) {
            cout << "          edge length inside polygon: " << edge_len_within_polygon << " total perim: " << P1 << " charge within poly " << Q1 << endl;
        }
        
        return Vec(
            F1_charge_within_poly * _E_x,
            F1_charge_within_poly * _E_y
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
        
        _E_y = E_y;
        _E_x = E_x;
        
        _poly_zone.set_points(polygon_vertices, verbose);
    }
    
    bool ForceEFieldOnCellBoundary::vec_outside_polygon_bbox(const Vec& v)
    {
        bool vec_inside_bounding_box = (
            v.x <= _poly_zone.xmax() &&
            v.x >= _poly_zone.xmin() &&
            v.y <= _poly_zone.ymax() &&
            v.y >= _poly_zone.ymin()
        );
        
        return ! vec_inside_bounding_box;
    }
    
    bool ForceEFieldOnCellBoundary::enabled_for_faceidx(int fid, bool verbose)
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
     
    Vec ForceEFieldOnCellBoundary::compute_he_force(const Vertex<Property>& vert, const HalfEdge<Property>& he, bool verbose)
    {
        
        if (verbose) {
            cout << "        ForceEFieldOnCellBoundary::compute_he_force - finding force by he " << he.idx() << endl;
        }
        
        Vec v_r = vert.data().r;
        Vec vto_r = he.to()->data().r;
        
        if (vec_outside_polygon_bbox(v_r) && vec_outside_polygon_bbox(vto_r)) {
            return Vec(0.0, 0.0);
        }
        
        double edge_len_within_polygon = _poly_zone.cell_edge_intersection(v_r, vto_r, verbose);
        
        const Face<Property>& f   = *(he.face());         // cell to the right of the half edge
        const Face<Property>& fp  = *(he.pair()->face()); // pair cell (opposite side of the same junction)
        
        // Total force on vertex by the half edge, x and y
        // double hev_force_x = 0.0; 
        Vec hev_force(0.0,0.0);
        
        if (enabled_for_faceidx(f.id, verbose)) {
            if (verbose) {
                cout << "          force enabled for fid " << f.id << ", calculating electric force" << endl;
            }
            
            hev_force += compute_he_force_for_face(vert, he, f, edge_len_within_polygon, verbose);
        } else {
            if (verbose) {
                cout << "          force not enabled for fid " << f.id << ", skipping..." << endl;
            }
        }
        if (enabled_for_faceidx(fp.id, verbose)) {
            if (verbose) {
                cout << "          force enabled for fid " << fp.id << ", calculating electric force" << endl;
            }
            
            hev_force += compute_he_force_for_face(vert, he, fp, edge_len_within_polygon, verbose);
        } else {
            if (verbose) {
                cout << "          force not enabled for fid " << fp.id << ", skipping..." << endl;
            }
        }
        
        if (verbose) {
            cout << "        ForceEFieldOnCellBoundary::compute_he_force - force on vertex " << vert.id << " by he " << he.idx() << " is :" << hev_force.x << ", " << hev_force.y << endl;
        }
        
        return hev_force;
    }
}
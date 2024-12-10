import math
import numpy as np
import json

class HexagonalCellMesh:
    
    @classmethod
    def get_cell_row_xpositions(cls, side_length, n_cells_in_row, row_is_odd):
        if (n_cells_in_row % 2 != 1 and row_is_odd) or (n_cells_in_row % 2 != 0 and not row_is_odd):
            raise Exception(
                "Odd rows should have odd number of cells, and even should have even number." +
                "But in this row, n_cells={}, oddness={}".format(n_cells_in_row, row_is_odd))
        
        if row_is_odd:
            middle_idx = (n_cells_in_row//2 + 1)
        else:
            middle_idx = n_cells_in_row/2 - 1/2
        
        # The cells are staggered like bricks, so they are spaced within a row by 2*side_length
        cells_in_row_spacing = 2*side_length
        
        xpositions = []
        for cell_idx in range(n_cells_in_row):
            cell_xpos = (cell_idx - middle_idx) * cells_in_row_spacing
            xpositions.append(cell_xpos)
        
        return xpositions
        
        
        
    
    # def generate_vertex_row(n_cells_horizontal)
    @classmethod
    def generate_vertex_row_untrimmed(cls, side_length, row_is_odd, num_cells_wide, verbose=False):
        # side
        if num_cells_wide % 2 != 1:
            raise ValueError("The tissue should be an odd number of cells wide.")
        
        vertex_x_positions = []
        
        # Two possibilities
        #########
        # Odd rows look like this - using parentheses to show middle position.:
        #    X------X--X--()--X--X------X
        #       0         1         2
        # Even rows look like this:
        #      X--X------X()X------X--X
        #       0         1         2
        
        #########
        # Odd rows look like this - using parentheses to show middle position.:
        #      X--X------X--X--()--X--X------X--X
        #       0         1         2         3
        # Even rows look like this:
        #    X------X--X------X()X------X--X------X
        #       0         1         2         3
        
        # Look at every other 'cell position', make two vertices for it.
        # range(0, num_cells_wide, 2)
        
        cell_indices_to_build_vtxs_from = list(range(0, num_cells_wide, 2))
        first_build_cell_is_leftmost = len(cell_indices_to_build_vtxs_from) % 2 == 1
        
        build_cell_vertices_wide = not (row_is_odd ^ first_build_cell_is_leftmost)
        
        vertex_x_positions = []
        
        # vertices_close_to_cells = not ( (len(cell_indices_to_build_vtxs_from) % 2 == 1) ^ row_is_odd)
        for cell_idx in cell_indices_to_build_vtxs_from:
            cell_xpos = (cell_idx - (num_cells_wide // 2 + 1)) * side_length * 1.5
            # print()
            if build_cell_vertices_wide:
                left_vtx_x = cell_xpos - side_length
                right_vtx_x = cell_xpos + side_length
            else:
                left_vtx_x = cell_xpos - side_length/2
                right_vtx_x = cell_xpos + side_length/2
            
            vertex_x_positions.append(left_vtx_x)
            vertex_x_positions.append(right_vtx_x)
        
        return vertex_x_positions
        
    
    @classmethod
    def build_tissue(cls, side_length, box_lx, box_ly, verbose=False):
        hex_h = side_length*math.sqrt(3)
        hex_w = side_length*2
        
        # The 'first' row has the height of one hexagon - each additional one adds half a hexagon height.
        # so...
        # tissue_height = (num_rows/2 + 1/2)*hexagon_height
        num_cell_rows = math.floor(2*((box_ly / hex_h) - 1/2))
        # Round down to an odd number - things will be easier if we have one cell at the origin, and build
        # outwards. otherwise, we might end up with asymettrical tissue and stuff
        num_cell_rows = ((num_cell_rows - 1)//2)*2 + 1
        
        
        # num_cell_columns = 
        #### Finding number of cell columns
        # The first cell has width hex_w. Each additional hexagon adds width 1.5*hex_w.
        # So for n_cells, total width is = hex_w + (n_cells - 1)*1.5*hex_w
        #                      width_tot = hex_w*((3/2)*n_cells - 1/2)
        #                      width_tot / (hex_w*3/2) +1/2 = n_cells
        # So finding n_cells:
        # print("Box lx={}   box_lx/(hex_w*(3/2))={}".format(box_lx, box_lx/(hex_w*3/2)))
        num_cell_columns = math.floor((box_lx / (hex_w*(3/2))) + (1/2))
        
        num_cell_columns = ( (num_cell_columns - 1)//2 )*2 + 1
        
        if verbose:
            print("hex_w={} hex_h={}".format(hex_w,hex_h))
            print("Number of cell columns: {}".format(num_cell_columns))
            print("Number of cell rows: {}".format(num_cell_rows))
        if num_cell_columns % 2 != 1:
            raise Exception("Should be odd number of columns... something is wrong")
        if num_cell_rows % 2 != 1:
            raise Exception("Should be odd number of rows... something is wrong.")
        # print("Number of cell columns: {}".format(num_cell_columns))
        
        # The middle row is at origin, and the others are spaced apart from it
        cell_row_ypositions = (hex_h/2)*(np.arange(num_cell_rows) - num_cell_rows//2)
        
        # print("Cell row y positions: {}".format(cell_row_ypositions))
        
        # vtx_row_positions = (hex_h/2)*(np.arange(num_cell_rows + 2) - (num_cell_rows + / // 2)
        vtx_rows = []
        for vtx_row_idx in range(num_cell_rows + 2):
            # The middle cell row is "odd" - has cell in it's center.
            # So, the row of vertices that are along the middle of that row are 'odd'.
            #  - this middle row of vertices has index  n_cell_rows//2 + 1
            # Since vertex rows alternate being 'even' or 'odd',
            # a row with index i is odd, if (i - (n_cell_rows//2 + 1)) % 2 == 0
            # print(vtx_row_idx)
            vtx_row_is_odd = (vtx_row_idx - (num_cell_rows//2 + 1)) % 2 == 0
            
            vtx_x_positions = cls.generate_vertex_row_untrimmed(
                side_length=side_length,
                row_is_odd=vtx_row_is_odd,
                num_cells_wide=num_cell_columns,
            )
            if verbose:
                print("# of vertex x positions generated: {}".format(len(vtx_x_positions)))
            vtx_row_y = (vtx_row_idx - (num_cell_rows/2)) * (hex_h/2)
            # vtx_y_positions = 
            # row_vertices
            vertices_in_row = [ {"x": vtx_x_positions[i], "y": vtx_row_y} for i in range(num_cell_columns + 1) ]
            
            
            vtx_rows.append(vertices_in_row)
        
        # There should be no duplicates
        # @TODO remove, slow
        # all_vertices_to_examine = 
        # for 
        
        ##### Calculate the number of cells in odd or even rows
        num_non_center_cell_columsn = num_cell_columns - 1
        n_cells_in_ODD_rows = 1 + (num_non_center_cell_columsn // 4)*2
        n_cells_in_EVEN_rows = num_cell_columns - n_cells_in_ODD_rows
        
        
        
        # cell_objs = {}
        cells_map = {}
        # vertices_map = {}
        # reverse_vertices_map = {}
        # vertex_objs = {}
        cell_id = 0
        vtx_id = 0
        
        ##############
        # Assign cells to their corresponding vertices
        if verbose:
            print("N cells in odd rows: {}".format(n_cells_in_ODD_rows))
            print("N cells in even rows: {}".format(n_cells_in_EVEN_rows))
            print("BUILDING CELL ROWS")
        for c_row_idx in range(num_cell_rows):
            # middle_cell
            # The middle row is odd, and there is odd number of rows -
            #  so if there are 5 rows, then 0,2,4 are odd, 1,3 are even.
            row_is_odd = abs(c_row_idx - num_cell_rows//2) % 2 == 0
            
            if row_is_odd:
                n_cells_in_row = n_cells_in_ODD_rows
            else:
                n_cells_in_row = n_cells_in_EVEN_rows
            
            # print("Tot number of cell columns: {}".format(num_cell_columns))
            # print()
            if verbose:
                print("   Row index={} row_is_odd={} n_cells_in_row={}".format(c_row_idx, row_is_odd, n_cells_in_row))
            
            cell_centers_x = cls.get_cell_row_xpositions(side_length, n_cells_in_row, row_is_odd)
            # cell_col
            center_col_idx = num_cell_columns // 2
            if row_is_odd:
                cell_column_indices = [i for i in range(num_cell_columns) if abs(i - center_col_idx) % 2 == 0]
                # cell_column_indices = list(np.arange(n_cells_in_row) * 2)
            else:
                cell_column_indices = [i for i in range(num_cell_columns) if abs(i - center_col_idx) % 2 != 0]
                # cell_column_indices = list(np.arange(n_cells_in_row) * 2 + 1)
            # print(cell_centers_x)
            
            
            # cells_map = {}
            # vertices_map = {}
            if len(cell_centers_x) != len(cell_column_indices):
                print("Cell center x values: {}".format(cell_centers_x))
                print("Cell column indices: {}".format(cell_column_indices))
                raise ValueError("Expected number of cell centers to match number of cell column indices")
            
            for cell_x_pos, c_col_idx in zip(cell_centers_x, cell_column_indices):
            # for cell_row_idx in 
                # for cell_y_pos
                try:
                    cell_vertices = [
                        vtx_rows[c_row_idx][c_col_idx],
                        vtx_rows[c_row_idx][c_col_idx + 1],
                        vtx_rows[c_row_idx + 1][c_col_idx + 1],
                        vtx_rows[c_row_idx + 2][c_col_idx + 1],
                        vtx_rows[c_row_idx + 2][c_col_idx],
                        vtx_rows[c_row_idx + 1][c_col_idx],
                    ]
                except:
                    print("Cell column indices: {}".format(cell_column_indices))
                    print("C_row_idx={} c_col_idx={}".format(c_row_idx, c_col_idx))
                    raise
                
                matching_vertex_ids = []
                for cell_vtx in cell_vertices:
                    # print(cell_vtx)
                    if "id" not in cell_vtx:
                        cell_vtx["id"] = vtx_id
                        cell_vtx["cells"] = []
                        vtx_id = vtx_id + 1
                    
                    cell_vtx["cells"].append(cell_id)
                    matching_vertex_ids.append(cell_vtx["id"])
                    
                # Add neighbors to the vertices - this represents edges
                for vi in range(len(cell_vertices)):
                    this_id = cell_vertices[vi]["id"]
                    prev_pos = vi - 1 if vi >= 1 else len(cell_vertices) - 1
                    next_pos = vi + 1 if vi + 1 < len(cell_vertices) else 0
                    prev_id = cell_vertices[prev_pos]["id"]
                    next_id = cell_vertices[next_pos]["id"]
                    
                    if prev_id == this_id or next_id == this_id:
                        raise Exception("??????")
                    
                    # cell_vertices["vi"]
                    if "neighbours" not in cell_vertices[vi]:
                        cell_vertices[vi]["neighbours"] = []
                    
                    if prev_id not in cell_vertices[vi]["neighbours"]:
                        cell_vertices[vi]["neighbours"].append(prev_id)
                    if next_id not in cell_vertices[vi]["neighbours"]:
                        cell_vertices[vi]["neighbours"].append(next_id)
                    # cell_vertices[vi]
                
                cells_map[cell_id] = {
                    "vertices": matching_vertex_ids,
                    "center_x": cell_x_pos,
                    "center_y": cell_row_ypositions[c_row_idx],
                }
                if verbose:
                    print("    Adding cell '{cell_id}'... column index={cell_colidx}".format(
                        cell_id=cell_id,
                        cell_colidx=c_col_idx,
                        ))
                
                cell_id += 1
        
        vertices_map = {}
        for vtx_row in vtx_rows:
            for vtx in vtx_row:
                if "id" in vtx:
                    vertices_map[vtx["id"]] = {
                        "x": vtx["x"],
                        "y": vtx["y"],
                        "cells": vtx["cells"],
                        "neighbours": vtx["neighbours"],
                    }
        
        # Remove unused vertices from the neighbor lists
        # Also decide which vertices are on boundary
        # @TODO - these should be set in the fucntion that defines the boundary! Not separately.
        for vtx_id, vtx in vertices_map.items():
            vtx["neighbours"] = [n_id for n_id in vtx["neighbours"] if n_id in vertices_map]
            
            # # vertices in the interior will be part of three faces
            # if len(vtx["cells"]) == 3:
            #     vtx["is_boundary"] = False
            # elif len(vtx["cells"]) < 3:
            #     vtx["is_boundary"] = True
            # else:
            #     print()
            #     raise Exception("Vertex should not be part of more than three cells")
        
        return cells_map, vertices_map
        
        
        
    def __init__(
        self,
        side_length,
        box_lx,
        box_ly,
    ):
        # raise NotImplementedError()
        
        # center_cell = 
        # cent
        # centself.
        self.side_length = side_length
        self.box_lx = box_lx
        self.box_ly = box_ly
        self.cells_map = None
        self.vertices_map = None
        
        self.build_cells()
    
    def build_cells(self,verbose=False):
        self.cells_map, self.vertices_map = self.build_tissue(self.side_length, self.box_lx, self.box_ly)    
        # raise Exception("Have not excuded the corners for even rows that are first or lat row")
        if verbose:
            print("BUILT CELLS")
            print(self.cells_map)
            print(self.vertices_map)
    
    def set_all_A0(self, A0_value):
        for cell_id in self.cells_map:
            self.cells_map[cell_id]["A0"] = A0_value

    def set_all_P0(self, P0_value):
        for cell_id in self.cells_map:
            self.cells_map[cell_id]["P0"] = P0_value
        
    @classmethod
    def _vtx_angle_to(cls, vtx1, vtx2):
        x1, y1 = vtx1["x"], vtx1["y"]
        x2, y2 = vtx2["x"], vtx2["y"]
        
        if x2-x1 == 0 and y1-y2 == 0:
            print("x1,y1:{},{}     x2,y2:{},{}".format(x1,y1,x2,y2))
            raise Exception("Coincident vertices")
        r1 = x1 + 1j*y1
        r2 = x2 + 1j*y2
        
        r_between = r2 - r1
        
        ang = np.angle(r_between)
        # numpy returns (-pi, pi], we want (0, 2*pi]
        if ang <= 0:
            ang += 2*np.pi
        
        # if it's a hair over now, set back to 2*np.pi
        if ang > 2*np.pi:
            ang = 2*np.pi
        
        return ang
        
    
    @classmethod
    def _angle_subtract(cls, angle1, angle2):
        """ Assumes angle1 and angle 2 are between 0 and 2*pi """
        diff_ang = angle1 - angle2
        
        ### Put both angles in 0 to 2 pi range
        while angle1 < 0:
            angle1 += 2*np.pi
        while angle1 > 2*np.pi:
            angle1 -= 2*np.pi
        while angle2 < 0:
            angle2 += 2*np.pi
        while angle2 > 2*np.pi:
            angle2 -= 2*np.pi
        # if angle1 < 0 or angle2 < 0 or angle1 > 2*np.pi or angle2 > 2*np.pi:
        #     raise ValueError("expected angles between 0 and 2pi. invalid args {}, {}".format(angle1, angle2))
        
        if angle2 > angle1:
            angle2 -= 2*np.pi
        # diff_angle = angle1 - angle2
        # raise NotImplementedError("This isn't doing what I think it's doing...")
        return angle1 - angle2
        
    
    @classmethod
    def generate_boundary_face_export_obj(cls, vertices_map, cell_id, verbose=False):
        # @NOTE this assuems there is only one boundary, around the perimeter
        
        # Start at the leftmost vertex - should be on the boundary
        # starting_vertex_id = list(vertices_map.keys())[0]
        starting_vertex_id = None
        for vtx_id, vtx in vertices_map.items():
            
            if starting_vertex_id is None or (vtx["x"] < vertices_map[starting_vertex_id]["x"]):
                starting_vertex_id = vtx_id
                
        if verbose:
            print("Starting vertex id: ", starting_vertex_id)
            print(vertices_map[starting_vertex_id])
        
        
        boundary_vertices = [starting_vertex_id]
        # We pretend that the vertex before the leftmost one was pointing directly to the left, so 
        # that it will find a vertex on the boundary when it scans for the next link int he chain. Otherwise it might start working towards the interior of the mesh.
        last_angle = None
        
        
        while True:
            current_vtx_id = boundary_vertices[-1]
            prev_vtx_id = boundary_vertices[-2] if len(boundary_vertices) > 1 else None
            current_vtx = vertices_map[current_vtx_id]
            
            next_candidates = current_vtx["neighbours"]
            next_candidates = [vid for vid in next_candidates if vid != prev_vtx_id]
            try:
                neighbor_angles = [cls._vtx_angle_to(current_vtx, vertices_map[neigh_id]) for neigh_id in next_candidates]
            except:
                print("All vertices: ")
                for vtx_id, vtx in vertices_map.items():
                    print("  {}: {}".format(vtx_id, vtx))
                print("Current vertex '{}': {}".format(current_vtx_id, current_vtx))
                print("Neighbor ids (candidates): {}".format(next_candidates))
                print("Neighbors: {}".format([vertices_map[vid] for vid in next_candidates]))
                raise
            # neighbor_angl
            if verbose:
                print("current vertex id: {}".format(current_vtx_id))
                print("        position: {},{}".format(current_vtx['x'], current_vtx['y']))
                print("        next candidates: {}".format(next_candidates))
            
            if last_angle is None: # the first edge to find
                neighbor_diffs_from_left_direction = [cls._angle_subtract(n_ang, np.pi) for n_ang in neighbor_angles]
                # next_vtx_id = next_candidates
                selected_cand_idx = np.argmax(neighbor_diffs_from_left_direction)
                if verbose:
                    print("    Neighbor diffs from pi: {}".format(neighbor_diffs_from_left_direction))
            else:
                neighbor_diffs_from_last_angle = [cls._angle_subtract(n_ang, last_angle + np.pi) for n_ang in neighbor_angles]
                selected_cand_idx = np.argmax(neighbor_diffs_from_last_angle)
                if verbose:
                    print("    neighbor diffs from previous angles: {}".format(neighbor_diffs_from_last_angle))
            next_vtx_id = next_candidates[selected_cand_idx]
            next_angle = neighbor_angles[selected_cand_idx]
            
            # if next_vtx_id == 5:
            #     print(neigh)
                # print(next_candidates)
            if next_vtx_id in boundary_vertices:
                if next_vtx_id == boundary_vertices[0]:
                    # Then we're circled back
                    break
                else:
                    # Then we've come back on the boundary without closing the loop - no good
                    raise Exception("Expected loop to close! Ran into boundary without reaching initial vertex first.")
            
            boundary_vertices.append(next_vtx_id)
            last_angle = next_angle
        
        
        
        
        
        return {
            "A0": 0,
            "P0": 0,
            "id": cell_id,
            "nsides": len(boundary_vertices),
            "vertices": boundary_vertices,
            "outer": True,
            "type": "passive",
        }
    
        
    def build_vm_mesh(self, out_json_fp, verbose=False):
        if self.cells_map is None or self.vertices_map is None:
            raise Exception("Cannot build vm mesh until cells have been built.")
            
            
        #### Build faces
        vm_faces = []
        for cell_id in self.cells_map:
            cell = self.cells_map[cell_id]
            cell_vertices = list(cell["vertices"])
            
            if "A0" not in cell:
                print(cell)
                raise Exception("Cannot build vm mesh without attribute 'A0' set. Cell id {} missing A0".format(cell_id))

            if "P0" not in cell:
                print(cell)
                raise Exception("Cannot build vm mesh without attribute 'P0' set. Cell id {} missing P0".format(cell_id))
            
            
            vm_faces.append({
                "A0": cell["A0"],
                "P0": cell["P0"],
                "id": cell_id,
                "nsides": len(cell_vertices),
                "vertices": cell_vertices,
                "outer": False,
                "type": "passive",
                # "vertices": cell_ve
            })
            
        boundary_face = self.generate_boundary_face_export_obj(self.vertices_map, cell_id=len(self.cells_map), verbose=verbose)
        vm_faces.append(boundary_face)
        # for v in boundary_face['vertices']:
        #     vm_
        
        
        #### Build vertices
        vm_vertices = []
        
        # raise NotImplementedError("Need to replace non indexed ids with real indexed ids!")
        
        for vertex_id in self.vertices_map:
            vtx = self.vertices_map[vertex_id]
            
            vtx_neighbors = list(vtx["neighbours"])
            
            vtx_type = "regular"
            vtx_x, vtx_y = vtx["x"], vtx["y"]
            
            vtx_is_boundary = vertex_id in boundary_face['vertices']
            
            num_faces = 0
            for face in vm_faces:
                if vertex_id in face["vertices"]:
                    num_faces += 1
            
            vm_vertices.append({
                "boundary": vtx_is_boundary,
                "constraint": "none",
                "coordination": num_faces,
                "erased": False,
                "id": vertex_id,
                "neighbours": vtx_neighbors,
                "r": [vtx_x, vtx_y],
                "type": vtx_type,
            })
            
        ##### Reindex faces
        # reindex_vm_faces = []
        for i, f in enumerate(vm_faces):
            f["id"] = i
            
        
        old_to_new_vid_map = {}
        ##### Reindex vertices
        for i, v in enumerate(vm_vertices):
            old_idx = v["id"]
            v["id"] = i
            old_to_new_vid_map[old_idx] = i
        
        for f in vm_faces:
            f["vertices"] = [old_to_new_vid_map[old_i] for old_i in f["vertices"]]
        
        
        
        # return res
        data = {
            "mesh": {
                "box": {
                    "lx": self.box_lx,
                    "ly": self.box_ly,
                    "periodic": False,
                },
                "faces": vm_faces,
                "vertices": vm_vertices,
            }
        }
        if verbose:
            print("Box:")
            print("  " + str(data["mesh"]["box"]))
            print("Faces:")
            for f in data["mesh"]["faces"]:
                print(f)
            print("Vertices:")
            for v in data["mesh"]["vertices"]:
                print(v)
            print("Old to new vertex id map: {}".format(old_to_new_vid_map))
            # print(json.dumps(data, indent=4))
        
        with open(out_json_fp, "w") as f:
            json.dump(data, f, indent=4)

class BoxSelector:
    @staticmethod
    def get_vertices_in_box(vm_state, x1,x2,y1,y2, verbose=False):
        assert x2 >= x1
        assert y2 >= y2
        
        if verbose:
            print("Sellecting vertices within xrange ({},{}) yrange ({},{})".format(
                x1,x2,y1,y2
            ))
        
        vertex_ids_inside = []
        
        for vid, vgeom in vm_state.current_tissue_state().geometry().vertices().items():
            vx = vgeom.x()
            vy = vgeom.y()
            
            if x1 <= vx and vx <= x2 and y1 <= vy and vy <= y2:
                vertex_ids_inside.append(vid)
        
        if verbose:
            print("SELECTED VERTICES:")
            print(vertex_ids_inside)
        return vertex_ids_inside
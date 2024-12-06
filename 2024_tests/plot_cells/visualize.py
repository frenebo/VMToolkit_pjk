
import pandas as pd
import plotly.graph_objects as go


# class TissStructure


class TissueVisualizer:
    @classmethod
    def tiss_struct_from_json_obj(cls, json_o):
        if "mesh" not in json_o:
            raise ValueError("Expected to find mesh in json object")
        if "vertices" not in json_o["mesh"]:
            raise ValueError("Expected to find vertices in json['mesh'] object")
        
        # vtx_id_map = {}
        vertices = {}
        for  v in json_o["mesh"]["vertices"]:
            if v["erased"]:
                raise NotImplementedError()
                
            vertices[v["id"]] = {
                "r": v["r"],
            }
        
        faces = {}
        for f in json_o["mesh"]["faces"]:
            if "erased" in f and f["erased"]:
                raise NotImplementedError()
            if f["outer"]:
                print("SKIPPING outer face - boundary")
                continue
            
            # If not outer, and not erased, we want to graph it
            
            faces[f["id"]] = {
                "vtx_ids": f["vertices"],
                "type": f["type"],
            }
        
        # return pd.DataFrame()
        return {
            "vertices": vertices,
            "faces": faces,
        }
        
        
    
    def __init__(self, init_json_o):
        self.structure = self.tiss_struct_from_json_obj(init_json_o)
        
        # self.update_cell_vertices(init_config)
    
    def add_frame(self, frame_idx, json_o):
        pass
        
        
    
    # def update_cell_vertices(self, json_conf):
        # self.structure
        
        # self.import_json_obj(init_config)

import pandas as pd
import plotly.graph_objects as go

import json

import os
# import plotly.graph_objects as go
import numpy as np


class TissueVisualizer:
    @classmethod
    def tiss_struct_from_json_obj(cls, json_o,verbose=False):
        if "mesh" not in json_o:
            raise ValueError("Expected to find mesh in json object")
        if "vertices" not in json_o["mesh"]:
            raise ValueError("Expected to find vertices in json['mesh'] object")
        
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
                if verbose:
                    print("SKIPPING outer face - boundary")
                continue
            
            # If not outer, and not erased, we want to graph it
            
            faces[f["id"]] = {
                "vtx_ids": f["vertices"],
                "type": f["type"],
            }
        
        return {
            "vertices": vertices,
            "faces": faces,
        }
        





def plotly_coords_list_from_checkpoint_obj(json_o_f):
    tiss_struct = TissueVisualizer.tiss_struct_from_json_obj(json_o_f)

    xx = []
    yy = []

    face_ids = sorted(list(tiss_struct["faces"].keys()))
    
    for f_id in face_ids:
        f_vids = list(tiss_struct["faces"][f_id]["vtx_ids"])

        # First index is repeated to make closed shape
        for vtx_id in f_vids + [f_vids[0]]:
            vtx_r = tiss_struct['vertices'][vtx_id]['r']
            xx.append(vtx_r[0])
            yy.append(vtx_r[1])
            
        xx.append(None)
        yy.append(None)
        
    return xx,yy
    
    
def build_framedata(checkpoint_objs):
    
    frame_coords = []
    for ckpt_obj in checkpoint_objs:
        xx,yy = plotly_coords_list_from_checkpoint_obj(ckpt_obj)
        frame_coords.append([xx,yy])

    frame_coords = np.array(frame_coords,dtype=float)
    
    

    xlims = (np.nanmin(frame_coords[:,0]), np.nanmax(frame_coords[:,0]))
    ylims = (np.nanmin(frame_coords[:,1]), np.nanmax(frame_coords[:,1]))
    
    return {
        "plotly_cell_coords": frame_coords,
        # "vtx_pos_range": {
        "xlims": xlims,
        "ylims": ylims,
        "n_frames": len(checkpoint_objs),
        # }
    }
    
    

def make_animated_tissue_plot(frame_data):
    fig = go.Figure(layout=go.Layout(
        xaxis=dict(
            range=frame_data["xlims"]
        ),
        yaxis=dict(
            range=frame_data["ylims"]
        ),
        width=600,
        height=600,
    ))

    fig.add_trace(
        go.Scatter(
            x=frame_data["plotly_cell_coords"][0][0],
            y=frame_data["plotly_cell_coords"][0][1],
            mode="lines", 
            marker_size=8,
            fill="toself",
            line=dict(color="red", width=2)
        ),
    )

    frames = [
        go.Frame(
            data=[
                {
                    "x": frame_data["plotly_cell_coords"][ckpt_idx][0],
                    "y": frame_data["plotly_cell_coords"][ckpt_idx][1],
                }
            ],
            name=ckpt_idx,
        ) for ckpt_idx in range(len(frame_data["plotly_cell_coords"]))
    ]


    sliders_dict = {
        "active": 0,
        "yanchor": "top",
        "xanchor": "left",
        "currentvalue": {
            "font": {"size": 20},
            "prefix": "Ckpt:",
            "visible": True,
            "xanchor": "right"
        },
        "transition": {"duration": 100, "easing": "cubic-in-out"},
        "pad": {"b": 10, "t": 50},
        "len": 0.9,
        "x": 0.1,
        "y": 0,
        "steps": []
    }
    for ckpt_idx in range(frame_data["n_frames"]):
        slider_step = {
            "args": [
                [ckpt_idx],
                {
                    "frame": {"duration": 100, "redraw": False},
                    "mode": "immediate",
                    "transition": {"duration": 100},
                }
            ],
            "label": "{}".format(ckpt_idx),
            "method": "animate"
        }
        sliders_dict["steps"].append(slider_step)




    fig.update_layout(
        sliders=[sliders_dict],
        updatemenus=[{
        "type": "buttons",
        "x": 1.075,
        "y":0.99,
        "buttons": [{
                "label":'Play',
                "method":'animate',    
                "args":[
                    None,
                    {
                        "frame":{"duration":100, "redraw":False}, 
                        "transition":{"duration":100},
                        # "fromcurrent":True,
                        "mode":'immediate'
                    },
                ]
        }]
    }])
    
    if frame_data["ylims"][1] - frame_data["ylims"][0] > frame_data["xlims"][1] - frame_data["xlims"][0]:
        fig.update_xaxes(
            scaleanchor="y",
            scaleratio=1,
            )
    else:
        fig.update_yaxes(
            scaleanchor="x",
            scaleratio=1,
            )
    
    fig.update(frames=frames)
    fig.show()



def make_plotly_visualizer(tiss_ckpt_fps):
    checkpoint_objs = []
    for ckpt_fp in tiss_ckpt_fps:
        with open(ckpt_fp, "r") as f:
            ckpt_obj = json.load(f)
        checkpoint_objs.append(ckpt_obj)
        
    
    frame_data = build_framedata(checkpoint_objs)
    
    make_animated_tissue_plot(frame_data)
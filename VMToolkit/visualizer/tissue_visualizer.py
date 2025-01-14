import pandas as pd
import plotly.graph_objects as go
import shapely

import json
import math

import os
import numpy as np


class TissueVisualizer:
    @classmethod
    def tiss_struct_from_json_obj(cls, json_o,verbose=False):
        
        # if "mesh" not in json_o:
        #     raise ValueError("Expected to find mesh in json object")
        # if "vertices" not in json_o["mesh"]:
        #     raise ValueError("Expected to find vertices in json['mesh'] object")
        
        vertices = {}
        for vid, v in json_o["current_tissue_state"]["geometry"]["vertices"].items():
            # if v["erased"]:
            #     raise NotImplementedError()
                
            vertices[vid] = {
                "r": [ v["x"], v["y"] ],
            }
        
        faces = {}
        for fid, f in json_o["topology"]["cells"].items():
            # if "erased" in f and f["erased"]:
            #     raise NotImplementedError()
            # if f["outer"]:
            #     if verbose:
            #         print("SKIPPING outer face - boundary")
            #     continue
            
            # If not outer, and not erased, we want to graph it
            
            faces[fid] = {
                "vtx_ids": f["vertex_ids"],
                # "type": f["type"],
            }
        
        return {
            "vertices": vertices,
            "faces": faces,
        }
        





def plotly_coords_list_from_checkpoint_obj(json_o_f, vertices_to_highlight=None):
    tiss_struct = TissueVisualizer.tiss_struct_from_json_obj(json_o_f)

    xx = []
    yy = []
    vtx_colors = []

    face_ids = sorted(list(tiss_struct["faces"].keys()))
    
    for f_id in face_ids:
        f_vids = list(tiss_struct["faces"][f_id]["vtx_ids"])

        # First index is repeated to make closed shape
        for vtx_id in f_vids + [f_vids[0]]:
            vtx_r = tiss_struct['vertices'][vtx_id]['r']
            xx.append(vtx_r[0])
            yy.append(vtx_r[1])
            if vertices_to_highlight is not None and vtx_id in vertices_to_highlight:
                vtx_colors.append("blue")
            else:
                vtx_colors.append("black")
            
        xx.append(None)
        yy.append(None)
        
    return {
        "xx":xx,
        "yy":yy,
        "colors": vtx_colors,
    }
    
    
def build_framedata(checkpoint_objs, vertices_to_highlight=None):
    
    frame_coords = []
    # vtx_colors = 
    frame_vtx_colors = []
    for ckpt_obj in checkpoint_objs:
        f_res= plotly_coords_list_from_checkpoint_obj(ckpt_obj, vertices_to_highlight=vertices_to_highlight)
        xx,yy,vtx_colors = f_res["xx"], f_res["yy"], f_res["colors"]
        frame_coords.append([xx,yy])
        frame_vtx_colors.append(vtx_colors)

    frame_coords = np.array(frame_coords,dtype=float)
    
    

    xlims = (np.nanmin(frame_coords[:,0]), np.nanmax(frame_coords[:,0]))
    ylims = (np.nanmin(frame_coords[:,1]), np.nanmax(frame_coords[:,1]))
    
    return {
        "plotly_cell_coords": frame_coords,
        "f_vtx_colors": frame_vtx_colors,
        # "vtx_pos_range": {
        "xlims": xlims,
        "ylims": ylims,
        "n_frames": len(checkpoint_objs),
        # }
    }
    
    

def make_animated_tissue_plot(frame_data, fields_data):
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
            # fill="toself",
            marker_color=frame_data["f_vtx_colors"][0],
            # line=dict(color="red", width=2)
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
        "transition": {"duration": 200, "easing": "cubic-in-out"},
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
                    "frame": {"duration": 200, "redraw": False},
                    "mode": "immediate",
                    "transition": {"duration": 200},
                }
            ],
            "label": "{}".format(ckpt_idx),
            "method": "animate"
        }
        sliders_dict["steps"].append(slider_step)
        
    # arrow_annotations = []
    
    # for x0,y0,x1,y1 in zip(x_end, y_end, x_start, y_start):
    #     arrow = go.layout.Annotation(dict(
    #                     x=x0,
    #                     y=y0,
    #                     xref="x", yref="y",
    #                     text="",
    #                     showarrow=True,
    #                     axref="x", ayref='y',
    #                     ax=x1,
    #                     ay=y1,
    #                     arrowhead=3,
    #                     arrowwidth=1.5,
    #                     arrowcolor='rgb(255,51,0)',)
    #                 )
    #     list_of_all_arrows.append(arrow)
    
    arrow_annotations = []
    field_colors = ["rgb(245, 84, 66)", "rgb(255, 5, 176)", "green", "purple"]
    if len(fields_data) > len(field_colors):
        raise Exception("Need to add more colorss")
    for field_i, field_dat in enumerate(fields_data):
        # print(field_dat)
        field_plot_col = field_colors[field_i]
        
        verts_x, verts_y = np.array(field_dat["vertices"],dtype=float).T
        verts_x = list(verts_x) + [verts_x[0]]
        verts_y = list(verts_y) + [verts_y[0]]
        

        fig.add_trace(
            go.Scatter(
                x=verts_x,
                y=verts_y,
                mode="lines", 
                marker_size=8,
                fill="toself",
                # marker_color="red",
                line={
                    "width": 2,
                    "dash": "dot",
                    "color": field_plot_col,
                }
            ),
        )
        
        arrow_annotations.append(
            go.layout.Annotation(dict(
                ax=field_dat["center_x"],
                ay=field_dat["center_y"],
                xref="x", yref="y",
                text="",
                showarrow=True,
                axref="x", ayref="y",
                x=field_dat["center_x"] + field_dat["E_n_x"]*10,
                y=field_dat["center_y"] + field_dat["E_n_y"]*10,
                arrowhead=3,
                arrowwidth=5,
                arrowcolor=field_plot_col,
            ))
        )
        # field_dat[""]

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
            }],
        }],
        annotations=arrow_annotations,
    )
    
    # fig.add_trace(
    #     go.Scatter(x=)
    # )
    
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
    # fig.show()
    
    return fig

def build_fields_data(vmstate_fp):
    with open(vmstate_fp, "r") as f:
        vmstate_obj = json.load(f)
    
    elec_fields = []
    for fid, fspec in vmstate_obj['current_tissue_state']['forces'].items():
        if fspec['type'] == "electric_cell_boundary_force":
            elec_fields.append(fspec['field_spec'])
    # print(elec_fields)
    
    fields_data = []
    
    for elec_f_spec in elec_fields:
        zone_center = shapely.centroid(shapely.Polygon(elec_f_spec["zone_bounds"]["polygon_vertices"]))
        E_x = elec_f_spec["E_x"]
        E_y = elec_f_spec["E_y"]
        
        if E_x == 0 and E_y == 0:
            continue
            
        E_amp = math.sqrt(E_x**2 + E_y**2)
        
        E_n_x = E_x / E_amp
        E_n_y = E_y / E_amp
        
        fields_data.append({
            "vertices": elec_f_spec["zone_bounds"]["polygon_vertices"],
            "center_x": zone_center.x,
            "center_y": zone_center.y,
            "E_n_x": E_n_x,
            "E_n_y": E_n_y,
        })
    
    return fields_data



def make_plotly_visualizer(tiss_ckpt_fps, init_vmstate_fp, vertices_to_highlight=None, janky_transitions=True):
    checkpoint_objs = []
    for ckpt_fp in tiss_ckpt_fps:
        with open(ckpt_fp, "r") as f:
            ckpt_obj = json.load(f)
        checkpoint_objs.append(ckpt_obj)
    
    if not janky_transitions:
        raise ValueError("Not having janky transitions is not supported")
        
    
    frame_data = build_framedata(checkpoint_objs, vertices_to_highlight=vertices_to_highlight)
    
    fields_data = build_fields_data(init_vmstate_fp)
    
    return make_animated_tissue_plot(
        frame_data,
        fields_data,
        # arrow_annotations,
    )
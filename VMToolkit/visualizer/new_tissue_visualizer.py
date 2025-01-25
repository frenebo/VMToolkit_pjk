
from ..sim.vm_state import (
    VMState,
    PixelatedElectricForceOnCellBoundary,
    UniformElectricForceOnCellBoundary,
)
import numpy as np
import plotly.graph_objects as go

def get_xy_lims_for_all_state_vertices(ckpt_states):
    all_x = []
    all_y = []
    for ck_state in ckpt_states:
        for vtxid, vtx_geo in ck_state.tissue_state().geometry().vertices().items():
            all_x.append(vtx_geo.x())
            all_y.append(vtx_geo.y())
    
    return (min(all_x), max(all_x)), (min(all_y), max(all_y))

def build_uniform_field_graph_objects(ck_state):
    for force_id in sorted(ck_state.tissue_state().forces().keys()):
        force = ck_state.tissue_state().forces()[force_id]
        
        if not isinstance(force, UniformElectricForceOnCellBoundary):
            continue
        
        raise NotImplementedError("Uniform field visualization not implemented yet")

def build_cell_mesh_graph_objects(ck_state):
    all_vtx_x = []
    all_vtx_y = []
    
    for cellidx, (cell_id, cell_top) in enumerate(sorted(ck_state.topology().cell_topologies().items())):
        if cellidx != 0:
            all_vtx_x.append(None)
            all_vtx_y.append(None)
        # cell_top = ck_state.topology().cell_topologies()[cell_id]
        
        vtx_ids = cell_top.vertex_ids()
        vertices_x = []
        vertices_y = []
        for vid in vtx_ids:
            vtx_geo = ck_state.tissue_state().geometry().vertices()[vid]
            vertices_x.append(vtx_geo.x())
            vertices_y.append(vtx_geo.y())
        all_vtx_x += vertices_x + [ vertices_x[0] ]
        all_vtx_y += vertices_y + [ vertices_y[0] ]
            
    return [
        go.Scatter(
            x=all_vtx_x,
            y=all_vtx_y,
            mode="lines", 
            marker_size=8,
            # fill="toself",
            marker_color='blue',
            # line=dict(color="red", width=2)
        )
    ]
        

def build_pixelated_field_graph_objects(ck_state, field_plot_col='green'):
    assert isinstance(ck_state, VMState)
    
    traces = []
    
    for force_id in sorted(ck_state.tissue_state().forces().keys()):
        force = ck_state.tissue_state().forces()[force_id]
        
        if not isinstance(force, PixelatedElectricForceOnCellBoundary):
            continue
        
        field_spec = force.field_spec()
        
        field_xmin = field_spec.grid_origin_x()
        field_xmax = field_spec.grid_origin_x() + field_spec.grid_spacing() * field_spec.grid_ncells_x()
        field_ymin = field_spec.grid_origin_y()
        field_ymax = field_spec.grid_origin_y() + field_spec.grid_spacing() * field_spec.grid_ncells_y()   
        
        corners = [
            ( field_xmin, field_ymin ),
            ( field_xmax, field_ymin ),
            ( field_xmax, field_ymax ),
            ( field_xmin, field_ymax ),
        ]
        
        plot_verts = np.array(corners + [ corners[0] ]) # close the loop
        traces.append(go.Scatter(
            x=plot_verts.T[0],
            y=plot_verts.T[1],
            mode="lines", 
            marker_size=8,
            fill="toself",
            # marker_color="red",
            line={
                "width": 2,
                "dash": "dot",
                "color": field_plot_col,
            }
        ))
        grid_lines_x = []
        grid_lines_y = []
        # vertical lines
        for col_idx in range(field_spec.grid_ncells_x()):
            col_x = field_xmin + col_idx * field_spec.grid_spacing()
            grid_lines_x += [
                col_x,
                col_x,
                None
            ]
            grid_lines_y += [
                field_ymin,
                field_ymax,
                None
            ]
        # horizontal lines
        for row_idx in range(field_spec.grid_ncells_y()):
            col_y = field_ymin + row_idx * field_spec.grid_spacing()
            grid_lines_x += [
                field_xmin,
                field_xmax,
                None
            ]
            grid_lines_y += [
                col_y,
                col_y,
                None
            ]
        # grid_lines_coords = np.array([grid_lines_x, grid_lines_y])
        # print(grid_lines_coords)
        traces.append(go.Scatter(
            x=grid_lines_x,
            y=grid_lines_y,
            mode='lines',
            # marker_color="red",
            line={
                "width": 1,
                "dash": "dot",
                "color": field_plot_col,
            }
        ))
            
    return traces

def build_force_arrow_graph_objects(ck_state, force_id):
    if ck_state.sim_current_stats() is None:
        return []
    
    vtx_current_forces = ck_state.sim_current_stats().forces_stats().current_vertex_forces().vertex_force_stats_by_forceid()[force_id]
    
    force_scale = 5
    
    arrows_x = []
    arrows_y = []
    for vtx_id, (force_x, force_y) in sorted(vtx_current_forces.force_experienced_by_vertex_id().items()):
        vtx_geo = ck_state.tissue_state().geometry().vertices()[vtx_id]
        
        arrows_x.append(vtx_geo.x())
        arrows_y.append(vtx_geo.y())
        arrows_x.append(vtx_geo.x() + force_x * force_scale)
        arrows_y.append(vtx_geo.y() + force_y * force_scale)
        arrows_x.append(None)
        arrows_y.append(None)
    
    return [go.Scatter(
        x=arrows_x,
        y=arrows_y,
        marker={
            "size": 10,
            "symbol": 'arrow-bar-up',
            "angleref": 'previous',
            "color": 'black',
        }
    )]
    

def build_data_for_frame(ck_state):
    pixelated_field_traces = build_pixelated_field_graph_objects(ck_state)
    
    mesh_traces = build_cell_mesh_graph_objects(ck_state)
    
    left_force_arrow_traces = build_force_arrow_graph_objects(ck_state, 'left_forcing_field')
    
    return mesh_traces + pixelated_field_traces + left_force_arrow_traces


def generate_frame_slider(n_frames):
    steps = []
    for ckpt_idx in range(n_frames):
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
        steps.append(slider_step)
        
    return {
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
        "steps": steps
    }
    
def create_plotly_figure(ckpt_states, verbose=False):
    xlims, ylims = get_xy_lims_for_all_state_vertices(ckpt_states)
    
    all_frames_data = [ build_data_for_frame(ck_state) for ck_state in ckpt_states ]
    all_frames = [go.Frame(data=frame_dat, name=frame_idx) for frame_idx, frame_dat in enumerate(all_frames_data)]
    
    print("length of all frames data: ", len(all_frames))
    # print(all_frames)
    
    fig = go.Figure(
        data = all_frames_data[0],
        layout=go.Layout(
            xaxis=dict(
                range=xlims,
            ),
            yaxis=dict(
                range=ylims,
            ),
            width=600,
            height=600,
            sliders=[
                generate_frame_slider(len(ckpt_states))
            ],
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
        ),
        frames=all_frames,
    )
    
    
    if ylims[1] - ylims[0] > xlims[1] - xlims[0]:
        fig.update_xaxes(
            scaleanchor="y",
            scaleratio=1,
            )
    else:
        fig.update_yaxes(
            scaleanchor="x",
            scaleratio=1,
            )
    
    return fig

def new_visualize_simulation(ckpt_states, verbose=False):
    return create_plotly_figure(ckpt_states, verbose=verbose)
    

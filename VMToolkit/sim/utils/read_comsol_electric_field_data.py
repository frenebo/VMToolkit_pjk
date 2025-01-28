

# import os
import numpy as np

from ..vm_state import PixelatedFieldSpec



def _read_comsol_pt_data(com_text):
    all_x = []
    all_y = []
    all_E = []
    
    for line in com_text.split("\n"):
        if line.startswith("%"):
            continue
        if len(line) == 0:
            continue
        xtext, ytext, E_x_text, E_y_text, E_z_text = line.split()
        # print(line.split())
        xcoord = float(xtext)
        ycoord = float(ytext)
        E_x = None if E_x_text == "NaN" else float(E_x_text)
        E_y = None if E_y_text == "NaN" else float(E_y_text)
        E_z = None if E_z_text == "NaN" else float(E_z_text)
        # print(E_x)
        if E_x is None or E_y is None or E_z is None:
            continue
        all_x.append(xcoord)
        all_y.append(ycoord)
        E = {
            "Ex": E_x,
            "Ey": E_y,
            "Ez": E_z,
        }
        all_E.append(E)

    return {
        "x_vals": all_x,
        "y_vals": all_y,
        "E_vals": all_E,
    }

def _calculate_grid_spacing(xvals, yvals):
    unique_sorted_xvals = np.sort(np.unique(xvals))
    unique_sorted_yvals = np.sort(np.unique(yvals))

    xspaces = unique_sorted_xvals[1:] - unique_sorted_xvals[:-1]
    yspaces = unique_sorted_yvals[1:] - unique_sorted_yvals[:-1]

    grid_spacing = xspaces[0]

    if not np.all(np.isclose(np.full_like(xspaces, grid_spacing), xspaces)):
        raise ValueError("X spaces don't all match: {}".format(xspaces))
    if not np.all(np.isclose(np.full_like(yspaces, grid_spacing), yspaces)):
        raise ValueError("Y and X spaces don't match: {} & {}".format(xspaces, yspaces))

    return grid_spacing


def read_comsol_electric_field_data(com_text, electric_field_force_multiplier=1.0, verbose=False):
    pointwise_field_vals = _read_comsol_pt_data(com_text)
    grid_spacing = _calculate_grid_spacing(
        pointwise_field_vals["x_vals"],
        pointwise_field_vals["y_vals"],
    )
    # print(grid_spacing)
    
    # the smallest x and y coordinates for which the field is defined
    xpoint_min = min(pointwise_field_vals["x_vals"])
    ypoint_min = min(pointwise_field_vals["y_vals"])
    
    
    xpoint_max = max(pointwise_field_vals["x_vals"])
    ypoint_max = max(pointwise_field_vals["y_vals"])
    
    n_pixels_x = round((xpoint_max - xpoint_min)/grid_spacing + 1)
    n_pixels_y = round((ypoint_max - ypoint_min)/grid_spacing + 1)
    
    field_data = [[[0,0] for col in range(n_pixels_y)] for row in range(n_pixels_x)] #np.zeros((n_pixels_x, n_pixels_y, 2), dtype=float)
    
    
    for x, y, E in zip(
        pointwise_field_vals["x_vals"],
        pointwise_field_vals["y_vals"],
        pointwise_field_vals["E_vals"],
        ):
        grid_x = round( (x - xpoint_min)/grid_spacing )
        grid_y = round( (y - ypoint_min)/grid_spacing )
        
        field_data[grid_x][grid_y][0] = E["Ex"] * electric_field_force_multiplier
        field_data[grid_x][grid_y][1] = E["Ey"] * electric_field_force_multiplier
        
    
    # the field values are defined at the CENTERS of the pixels - so the pixel origins are shifted by half a pixel in each direction
    grid_xorigin = xpoint_min - grid_spacing/2
    grid_yorigin = ypoint_min - grid_spacing/2
    
    if verbose:
        print("Returning a pixelated field spec with origin {},{}".format(grid_xorigin, grid_yorigin))
        print("   grid spacing={}".format(grid_spacing))
        print("   npixels x,y={},{}".format(n_pixels_x, n_pixels_y))
    
    return PixelatedFieldSpec(
        grid_origin_x=grid_xorigin,
        grid_origin_y=grid_yorigin,
        grid_spacing=grid_spacing,
        grid_ncells_x=n_pixels_x,
        grid_ncells_y=n_pixels_y,
        field_data=field_data,
    )

    # # return ComsolEFieldData(
    # )
    

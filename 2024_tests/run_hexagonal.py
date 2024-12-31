import numpy as np
import argparse
import json
import os
from line_profiler import profile

from simcode.theoretical_model.find_hexagon_rest_area import HexagonalModel, find_hexagon_rest_area
from VMToolkit.config_builder.open.honeycomb_lattice import HoneycombLattice

from VMToolkit.VMAnalysis.utils.HalfEdge import Mesh

from simcode.tissue_builder.hexagonal import HexagonalCellMesh
from simcode.sim_model.sim_model import SimModel
from simcode.sim_model.vm_state import (
    VMState, SimulationSettings, IntegratorSettings, 
    CellAreaForce, CellPerimeterForce, ConstantVertexForce,  CellGroup, VertexGroup,
    ElectricForceOnCellBoundary, EFieldSpecConstantPolygonRegion, PolygonSpec,
    
)
from simcode.sim_model.box_selector import BoxSelector


def make_forcing_field_rectangular(
    xmin,
    xmax,
    ymin,
    ymax,
    field_x,
    field_y,
    ):
    return ElectricForceOnCellBoundary(electric_field_spec=(
        EFieldSpecConstantPolygonRegion(
            E_x=field_x,
            E_y=field_y,
            zone_bounds=PolygonSpec(
                polygon_vertices=[
                    (xmin,ymin),
                    (xmax,ymin),
                    (xmax,ymax),
                    (xmin,ymax),
                ]
            )
        )
    ))

@profile
def do_stuff():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fpull", default=-0.2,type=float, help="Force to squeeze/stretch")
    parser.add_argument("--box_lx", default=30.0,type=float, help="X size of tissue")
    parser.add_argument("--box_ly", default=30.0,type=float, help="Y size of tissue")
    parser.add_argument("--field_strength", default=0.00, type=float, help="field strength")
    
    args = parser.parse_args()
    
    A0_model = 5
    P0_model = 6.0
    gamma = 0.1
    kappa = 0.1  
    
    analytical_predictions = HexagonalModel().find_rest_size_of_hexagon(
        A_0=A0_model,
        P_0=P0_model,
        K=kappa,
        Gamma=gamma,
    )
    
    rest_side_length = analytical_predictions["rest_side_length"]
    theoretical_rest_width = rest_side_length*2
    theoretical_rest_height = rest_side_length*2*np.sqrt(3)/2
    
    print("Theoretical rest width, height: {}, {}".format(theoretical_rest_width, theoretical_rest_height))
    
    
    cm = HexagonalCellMesh(
        side_length=rest_side_length,
        box_lx=args.box_lx,
        box_ly=args.box_ly,
        # verbose=True,
    )
    tiss_topology, tiss_init_state = cm.build_vm_state()#verbose=True)
    
    tiss_init_state.forces()["perim_f_all"] = CellPerimeterForce(gamma=gamma, lam=P0_model*gamma)
    tiss_init_state.forces()["area_f_all"] = CellAreaForce(A0=A0_model, kappa=kappa)
    
    field_size_multiplier = 0.2
    tiss_init_state.forces()["left_forcing_field"] = make_forcing_field_rectangular(
        xmin=2*field_size_multiplier*(-args.box_lx),
        xmax=2*field_size_multiplier*(-args.box_lx)/3,
        ymin=field_size_multiplier*(-args.box_ly),
        ymax=field_size_multiplier*(args.box_ly),
        field_x=0.000,
        field_y=args.field_strength,
    )
    tiss_init_state.forces()["right_forcing_field"] = make_forcing_field_rectangular(
        xmin=2*field_size_multiplier*(args.box_lx)/3,
        xmax=2*field_size_multiplier*(args.box_lx),
        ymin=field_size_multiplier*(-args.box_ly),
        ymax=field_size_multiplier*(args.box_ly),
        field_x=0.000,
        field_y=-args.field_strength,
    )
    
    tiss_init_state.cell_groups()["all"].force_ids().extend([
        "perim_f_all",
        "area_f_all",
        "left_forcing_field",
        "right_forcing_field",
    ])
    
    vm_initial_state = VMState(
        tiss_topology=tiss_topology,
        current_state=tiss_init_state,
        sim_settings=SimulationSettings(
            integrator_settings=IntegratorSettings(
                vertex_friction_gamma=0.1,
                step_dt=0.08,
            ),
        )
    )
    
    
    print("MODEL PARAMS")
    print("A0={A0}  P0={P0}  kappa={kappa} gamma={gamma}".format(A0=A0_model,P0=P0_model,kappa=kappa,gamma=gamma))
    print(analytical_predictions)
    print("theoretical rest width={}, height={}".format(theoretical_rest_width, theoretical_rest_height))


    # # #################################################################
    # #
    # # Create the initial configuration and read it
    # #
    # # #################################################################
    print("Instantiating SimModel")
    sim_model = SimModel(verbose=False)
    print("Running sim_model.load_from_json_state")
    with open("scratch/initial_vm_state.json", "w") as f:
        f.write(json.dumps(vm_initial_state.to_json()))
    sim_model.load_from_json_state(vm_initial_state.to_json())#, verbose=True)
    print("Finished sim_mode.load_from_json_state")
    
    checkpoint_fps = []
    
    ckpt_dir = "scratch"
    if not os.path.exists(ckpt_dir):
        raise Exception("could not find {}".format(ckpt_dir))
    
    for fn in os.listdir(ckpt_dir):
        if (fn.startswith("res") or fn.startswith("vmst_")) and fn.endswith(".json"):
            os.remove(os.path.join(ckpt_dir, fn))
    
    step_size = 1000     # Step counter in terms of time units
    N_checkpoints = 50
    
    for i in range(N_checkpoints):
        ckpt_fp = "scratch/res{}.json".format(str(i).zfill(3))
        vmstate_json = sim_model.vm_state_json()
        with open(ckpt_fp, "w") as f:
            f.write(json.dumps(vmstate_json))
        
        checkpoint_fps.append(ckpt_fp)
        sim_model.run_steps(step_size, do_time_force_computation=True)


if __name__ == "__main__":
    do_stuff()

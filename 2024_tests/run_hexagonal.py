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
    TopologySettings, T1TransitionSettings,
    
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
    
    parser.add_argument("--box_nsides_x", default=40,type=float, help="X size of tissue")
    parser.add_argument("--box_nsides_y", default=40,type=float, help="Y size of tissue")
    parser.add_argument("--field_strength", default=1.6, type=float, help="field strength")
    
    args = parser.parse_args()
    
    A0_model = 9
    P0_model = 5
    gamma = 1.0
    kappa = 1.0
    
    analytical_predictions = HexagonalModel().find_rest_size_of_hexagon(
        A_0=A0_model,
        P_0=P0_model,
        K=kappa,
        Gamma=gamma,
    )
    
    rest_side_length = analytical_predictions["rest_side_length"]
    theoretical_rest_width = rest_side_length*2
    theoretical_rest_height = rest_side_length*2*np.sqrt(3)/2
    print(analytical_predictions)
    print("Theoretical rest width, height: {}, {}".format(theoretical_rest_width, theoretical_rest_height))
    
    box_lx = args.box_nsides_x*rest_side_length
    box_ly = args.box_nsides_y*rest_side_length
    cm = HexagonalCellMesh(
        side_length=rest_side_length,
        box_lx=box_lx,
        box_ly=box_ly,
    )
    tiss_topology, tiss_init_state = cm.build_vm_state()#verbose=True)
    
    tiss_init_state.forces()["perim_f_all"] = CellPerimeterForce(gamma=gamma, lam=P0_model*gamma)
    tiss_init_state.forces()["area_f_all"] = CellAreaForce(A0=A0_model, kappa=kappa)
    
    field_size_multiplier = 0.08
    tiss_init_state.forces()["left_forcing_field"] = make_forcing_field_rectangular(
        xmin=2*field_size_multiplier*(-box_lx),
        xmax=2*field_size_multiplier*(-box_lx)*0.0,
        ymin=2*field_size_multiplier*(-box_ly),
        ymax=2*field_size_multiplier*(box_ly),
        # field_x=args.field_strength,
        field_x=0,
        field_y=args.field_strength,
        # field_y=args.field_strength,
    )
    tiss_init_state.forces()["right_forcing_field"] = make_forcing_field_rectangular(
        xmin=2*field_size_multiplier*(box_lx)*0.0,
        xmax=2*field_size_multiplier*(box_lx),
        ymin=2*field_size_multiplier*(-box_ly),
        ymax=2*field_size_multiplier*(box_ly),
        field_x=0,
        field_y=-args.field_strength,
    )
    # tiss_init_state.forces()["top_forcing_field"] = make_forcing_field_rectangular(
    #     ymin=2*field_size_multiplier*(-box_lx),
    #     ymax=2*field_size_multiplier*(-box_lx)*0.1,
    #     xmin=2*field_size_multiplier*(-box_ly),
    #     xmax=2*field_size_multiplier*(box_ly),
    #     field_x=0,
    #     field_y=-args.field_strength*2.0,
    #     # field_y=args.field_strength,
    # )
    # tiss_init_state.forces()["bottom_forcing_field"] = make_forcing_field_rectangular(
    #     ymin=2*field_size_multiplier*(box_lx)*0.1,
    #     ymax=2*field_size_multiplier*(box_lx),
    #     xmin=2*field_size_multiplier*(-box_ly),
    #     xmax=2*field_size_multiplier*(box_ly),
    #     field_x=0,
    #     field_y=args.field_strength*2.0,
    # )
    
    tiss_init_state.cell_groups()["all"].force_ids().extend([
        "perim_f_all",
        "area_f_all",
        "left_forcing_field",
        "right_forcing_field",
        # "top_forcing_field",
        # "bottom_forcing_field",
    ])
    
    vm_initial_state = VMState(
        tiss_topology=tiss_topology,
        current_state=tiss_init_state,
        sim_settings=SimulationSettings(
            integrator_settings=IntegratorSettings(
                vertex_friction_gamma=0.1,
                step_dt=0.012,
            ),
            topology_settings=TopologySettings(
                T1_transition_settings=T1TransitionSettings(
                    enabled=True,
                    min_edge_len=rest_side_length*0.2,
                    new_edge_len=rest_side_length*0.21,
                    # min_edge_len=0.2,
                    # new_edge_len=0.22,
                )
            )
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
    
    step_size = 100     # Step counter in terms of time units
    N_checkpoints = 100
    
    for i in range(N_checkpoints):
        ckpt_fp = "scratch/res{}.json".format(str(i).zfill(3))
        vmstate_json = sim_model.vm_state_json()
        with open(ckpt_fp, "w") as f:
            f.write(json.dumps(vmstate_json))
        
        checkpoint_fps.append(ckpt_fp)
        sim_model.run_steps(step_size, do_time_force_computation=True)
        if sim_model._check_topology_changed():
            print("TOP CHANGED _ changing step size")
            step_size = 1
            # break
    

if __name__ == "__main__":
    do_stuff()

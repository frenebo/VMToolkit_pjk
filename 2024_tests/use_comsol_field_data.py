import argparse
import math
import os
import numpy as np
import json

from VMToolkit.sim.utils import read_comsol_electric_field_data
from VMToolkit.theoretical_models import TheoreticalRegularHexModel
from VMToolkit.sim.tissue_builder import HexagonalCellMeshBuilder
from VMToolkit.sim import SimulationModel
from VMToolkit.sim.vm_state import (
    CellPerimeterForce,
    CellAreaForce,
    VMState,
    SimulationSettings,
    IntegratorSettingsDormandPrinceRungeKutta,
    TopologySettings,
    T1TransitionSettings,
    PixElecForceCellParam,
    PixelatedElectricForceOnCellBoundary,
)


def generate_vertex_model_energy_parameters(shape_param, Y_param, cell_equilibrium_area):
    init_A0 = 1
    init_P0 = shape_param * math.sqrt(init_A0)
    
    # Y = gamma_reg_value / (K_reg_value*A_r_reg_value)
    init_kappa = 1
    init_gamma = Y_param * init_kappa * init_A0
        
    init_elastic_props = TheoreticalRegularHexModel().find_elastic_props_of_hexagon(
        A_0_num=init_A0,
        P_0_num=init_P0,
        K_num=init_kappa,
        gamma_num=init_gamma,
        # verbose=True,
    )
    
    # print(init_elastic_props)
    
    # Correcting A0 and P0 to give the desired cell_equilibrium_area
    # area_correction_factor = cell_equilibrium_area / init_elastic_props["rest_area"]
    # print("Desired rest area: {}".format(cell_equilibrium_area))
    # print("Actual rest area: {}".format(init_elastic_props["rest_area"]))
    
    length_scale_correction_factor = math.sqrt(cell_equilibrium_area / init_elastic_props["rest_area"])
    # print("Length scale correction factor: {}".format(length_scale_correction_factor))
    # print("Area correction factor: {}".format(area_correction_factor))
    # step2_A0 = init_A0 * area_correction_factor
    # step2_P0 = init_P0 * np.sqrt(area_correction_factor)
    
    # unscaled_rest_area = ( (3 * (3**0.5))/2 ) * (init_elastic_props["rest_side_length"]**2)
    # print("Recalculated actual rest area: {}".format(unscaled_rest_area))    
    
    area_correction_factor = cell_equilibrium_area / init_elastic_props["rest_area"]
    step1_A0  = init_A0 * (length_scale_correction_factor ** 2)
    step1_P0 =  init_P0 * (length_scale_correction_factor)
    # print("SHAPE PARAM", step2_P0 / np.sqrt(step2_A0))
    
    # gamma units - Energy / length^2
    # kappa units - energy / length^4
    step1_gamma = init_gamma / (length_scale_correction_factor ** 2)
    step1_kappa = init_kappa / (length_scale_correction_factor ** 4)
    
    step1_elastic_props = TheoreticalRegularHexModel().find_elastic_props_of_hexagon(
        A_0_num=step1_A0,
        P_0_num=step1_P0,
        K_num=step1_kappa,
        gamma_num=step1_gamma,
    )
    
    if not np.isclose(step1_gamma / (step1_kappa * step1_A0), Y_param):
        raise ValueError("Y param should remain unchanged after change of length scale")
    
    return {
        "A0": step1_A0,
        "P0": step1_P0,
        "gamma": step1_gamma,
        "kappa": step1_kappa,
    }
    
def rest_side_length_for_rest_area(rest_area):
    
    # rest_area = (3*(math.sqrt(3))/2) * (rest_side_length ** 2)
    return math.sqrt(rest_area) / math.sqrt((3*(math.sqrt(3))/2))

def run_sim(args):
    cell_shape_param = args.shape_param
    cell_y_param = args.y_param
    cell_equilibrium_area = args.cell_equilibrium_area
    tissue_width = args.tissue_width
    tissue_height = args.tissue_height
    vertex_friction_gamma = args.vertex_friction_gamma
    init_integrator_dt = args.init_integrator_dt
    
    relative_dormand_prince_max_error = args.dormand_prince_max_error_rel
    relative_t1_min_edgelen = args.t1_min_edge_length_rel
    relative_t1_new_edgelen = args.t1_new_edge_length_rel
    
    electric_field_force_multiplier = args.electric_field_force_multiplier
    
    n_checkpoints = args.n_checkpoints
    ckpt_period = args.checkpoint_period
    
    outdir = args.outdir
    
    verbose = args.verbose
    # parser.add_argument("--dormand_prince_max_error_rel", type=float, default=0.002, help="max error allowed in the dormand prince runge kutta (relative to rest cell edge length)")
    # parser.add_argument("--t1_min_edge_length_rel", type=float, default=0.02, help="edge length at which t1 transition happens (relative to rest cell edge length)")
    # parser.add_argument("--t1_new_edge_length_rel", type=float, default=0.03, help="new edge length after t1 transition (relative to rest cell edge length)")
    
    # parser.add_argument("--init_integrator_dt", type=float, default=0.01, help="default dt for the runge kutta adaptive timestep")
    # parser.add_argument("--dormandprince_max_error", type=float, default=0.01, help="max error allowed in the dormand prince runge kutta")

    
    with open(args.comsoldata, "r") as f:
        electric_field_spec = read_comsol_electric_field_data(
            f.read(),
            electric_field_force_multiplier=electric_field_force_multiplier,
            verbose=True)
    
    rest_side_length = rest_side_length_for_rest_area(cell_equilibrium_area)
    
    dormand_prince_max_error = args.dormand_prince_max_error_rel * rest_side_length
    t1_min_edgelen = args.t1_min_edge_length_rel * rest_side_length
    t1_new_edgelen = args.t1_new_edge_length_rel * rest_side_length
    
    cm = HexagonalCellMeshBuilder(
        side_length=rest_side_length,
        box_lx=tissue_width,
        box_ly=tissue_height,
    )
    tiss_topology, tiss_init_state = cm.build_vm_state()
    
    energy_params = generate_vertex_model_energy_parameters(
        shape_param=cell_shape_param,
        Y_param=cell_y_param,
        cell_equilibrium_area=cell_equilibrium_area,
    )
    
    tiss_init_state.forces()["perim_f_all"] = CellPerimeterForce(
        gamma=energy_params["gamma"],
        lam=energy_params["P0"]*energy_params["gamma"],
    )
    tiss_init_state.forces()["area_f_all"] = CellAreaForce(
        A0=energy_params["A0"],
        kappa=energy_params["kappa"],
    )
    
    pix_efield_cell_params = {}
    for cell_id in tiss_init_state.cell_groups()["normal"].cell_ids():
        pix_efield_cell_params[cell_id] = PixElecForceCellParam(charge=1.0)
    
    tiss_init_state.forces()["comsol_field_force"] = PixelatedElectricForceOnCellBoundary(
        field_spec=electric_field_spec,
        cell_params=pix_efield_cell_params,
    )
    
    tiss_init_state.cell_groups()["normal"].force_ids().extend([
        "perim_f_all",
        "area_f_all",
        "comsol_field_force"
    ])
    
    vm_initial_state = VMState(
        tiss_topology=tiss_topology,
        current_tissue_state=tiss_init_state,
        sim_settings=SimulationSettings(
            integrator_settings=IntegratorSettingsDormandPrinceRungeKutta(
                vertex_friction_gamma=vertex_friction_gamma,
                init_dt=              init_integrator_dt,
                displacement_error_max=dormand_prince_max_error,
            ),
            topology_settings=TopologySettings(
                T1_transition_settings=T1TransitionSettings(
                    enabled=True,
                    min_edge_len=rest_side_length*relative_t1_min_edgelen,
                    new_edge_len=rest_side_length*relative_t1_new_edgelen,
                )
            )
        ),
        sim_current_stats=None,
    )
    
    for fn in os.listdir(outdir):
        if (fn.startswith("res") or fn.startswith("vmst_")) and fn.endswith(".json"):
            os.remove(os.path.join(outdir, fn))
            
            
    print("Instantiating SimulationModel")
    sim_model = SimulationModel(verbose=verbose)
    print("Running sim_model.load_from_json_state")
    
        
    with open(os.path.join(outdir, "initial_vm_state.json"), "w") as f:
        f.write(json.dumps(vm_initial_state.to_json(), indent=4))
    sim_model.load_from_json_state(vm_initial_state.to_json(), verbose=verbose)#, verbose=True)
    print("Finished sim_mode.load_from_json_state")
    
    
    
    
    ckpt_strnum_nchars = len(str(n_checkpoints - 1))
    
    tot_time_elaspsed = 0
    sim_model.update_vm_state_from_cpp_vm(verbose=verbose) # So we get what the instantaneous forces would be right at start
    for i in range(n_checkpoints):
        ckpt_filename = "res{}.json".format(str(i).zfill(ckpt_strnum_nchars))
        ckpt_fp = os.path.join(outdir, ckpt_filename)
        vmstate_json = sim_model.vm_state_json()
        # update_vm_state_from_cpp_vm
        with open(ckpt_fp, "w") as f:
            f.write(json.dumps(vmstate_json, indent=4))
        
        # checkpoint_fps.append(ckpt_fp)
        if verbose:
            print("Running with adaptive timestep")
        sim_model.run_with_adaptive_tstep(ckpt_period, do_time_force_computation=True, verbose=verbose)
        tot_time_elaspsed += ckpt_period
        print("Time elapsed: {}".format(tot_time_elaspsed))
        if sim_model._check_topology_changed():
            print("TOPOLOGY CHANGED")
    
def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument("--comsoldata", required=True, help="Path to a txt file with the field x,y,z values")
    parser.add_argument("--shape_param", type=float, required=True, help="Desired shape parameter - P0/sqrt(A0)")
    
    parser.add_argument("--y_param",type=float, required=True, help="Desired Y param - Gamma/(Kappa * A0)")
    parser.add_argument("--cell_equilibrium_area", type=float, default=16*0.00059180327, help="Desired equilibrium area for the cells to be, as determined by vertex model energy parameters")
    
    parser.add_argument("--tissue_width", type=float, default=6.0, help="Width of tissue in mm (6mm for conv)")
    parser.add_argument("--tissue_height", type=float, default=1.5, help="Height of tissue in mm (1.5mm for conv)")
    parser.add_argument("--vertex_friction_gamma", type=float, default=1.0, help="Vertex friction experienced")
    
    parser.add_argument("--init_integrator_dt", type=float, default=0.01, help="default dt for the runge kutta adaptive timestep")
    parser.add_argument("--dormand_prince_max_error_rel", type=float, default=0.002, help="max error allowed in the dormand prince runge kutta (relative to rest cell edge length)")
    parser.add_argument("--t1_min_edge_length_rel", type=float, default=0.02, help="edge length at which t1 transition happens (relative to rest cell edge length)")
    parser.add_argument("--t1_new_edge_length_rel", type=float, default=0.03, help="new edge length after t1 transition (relative to rest cell edge length)")
    
    parser.add_argument("--checkpoint_period", type=float, default=0.04, help="checkpoint period timelength")
    parser.add_argument("--n_checkpoints", type=int, default=20, help="How many checkpoints to do")
    
    parser.add_argument("--electric_field_force_multiplier", type=float, default=0.001, help="How much to scale up or down the electric field from comsol")
    
    parser.add_argument("--outdir", required=True)
    
    parser.add_argument("--verbose", action='store_true')
                # init_dt=              init_integrator_dt,
                # displacement_error_max=dormand_prince_max_error,
    # 61 cells per 0.19mm*0.19mm ->0.0361mm^2 -> cell area is 0.00059180327 mm^2
    
    args = parser.parse_args()
    
    run_sim(args)

if __name__ == "__main__":
    main()
    
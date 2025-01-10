import numpy as np
import argparse
import json
import math
import os
from dataclasses import dataclass

from simcode.theoretical_model.regular_hexagon_sympy_model import RegHexagonalModel, find_regular_hexagonal_rests_area
from simcode.tissue_builder.hexagonal import HexagonalCellMesh
from simcode.sim_model.sim_model import SimModel
from simcode.sim_model.vm_state import (
    VMState, SimulationSettings, IntegratorSettings, 
    CellAreaForce, CellPerimeterForce, ConstantVertexForce,  CellGroup, VertexGroup,
    ElectricForceOnCellBoundary, EFieldSpecConstantPolygonRegion, PolygonSpec,
    TopologySettings, T1TransitionSettings,
    IntegratorSettingsDormandPrinceRungeKutta,
)


def make_forcing_field_rectangular(
    xmin,
    xmax,
    ymin,
    ymax,
    field_x,
    field_y,
    ):
    return ElectricForceOnCellBoundary(
        electric_field_spec=EFieldSpecConstantPolygonRegion(
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
    )

def run_hexagonal(args, outdir):
    A0_model = args.param_A0
    P0_model = args.param_P0
    
    gamma = args.param_gamma
    kappa = args.param_kappa
    
    analytical_predictions = RegHexagonalModel().find_rest_size_of_hexagon(
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
    
    rest_area = (3*(math.sqrt(3))/2) * rest_side_length*rest_side_length
    
    box_lx = args.box_size_rel_x*math.sqrt(rest_area)
    box_ly = args.box_size_rel_y*math.sqrt(rest_area)
    
    cm = HexagonalCellMesh(
        side_length=rest_side_length,
        box_lx=box_lx,
        box_ly=box_ly,
    )
    tiss_topology, tiss_init_state = cm.build_vm_state()
    
    tiss_init_state.forces()["perim_f_all"] = CellPerimeterForce(gamma=gamma, lam=P0_model*gamma)
    tiss_init_state.forces()["area_f_all"] = CellAreaForce(A0=A0_model, kappa=kappa)
    
    field_size_multiplier = args.field_size_multiplier
    
    tiss_init_state.forces()["left_forcing_field"] = make_forcing_field_rectangular(
        xmin=(1/2)*field_size_multiplier*(-box_lx),
        xmax=(1/2)*field_size_multiplier*(-box_lx)*0.0,
        ymin=(1/2)*field_size_multiplier*(-box_ly),
        ymax=(1/2)*field_size_multiplier*(box_ly),
        field_x=args.field_strength,
        field_y=0,
    )
    tiss_init_state.forces()["right_forcing_field"] = make_forcing_field_rectangular(
        xmin=(1/2)*field_size_multiplier*(box_lx)*0.0,
        xmax=(1/2)*field_size_multiplier*(box_lx),
        ymin=(1/2)*field_size_multiplier*(-box_ly),
        ymax=(1/2)*field_size_multiplier*(box_ly),
        field_x=-args.field_strength,
        field_y=0,
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
    
    tiss_init_state.cell_groups()["normal"].force_ids().extend([
        "perim_f_all",
        "area_f_all",
        "left_forcing_field",
        "right_forcing_field",
        # "top_forcing_field",
        # "bottom_forcing_field",
    ])
    
    vm_initial_state = VMState(
        tiss_topology=tiss_topology,
        current_tissue_state=tiss_init_state,
        sim_settings=SimulationSettings(
            integrator_settings=IntegratorSettingsDormandPrinceRungeKutta(
                vertex_friction_gamma=args.vertex_friction_gamma,
                init_dt=args.init_integrator_dt,
                displacement_error_max=args.dormand_prince_max_error,
            ),
            topology_settings=TopologySettings(
                T1_transition_settings=T1TransitionSettings(
                    enabled=True,
                    min_edge_len=rest_side_length*args.t1_min_edge_len_rel,
                    new_edge_len=rest_side_length*args.t1_new_edge_len_rel,
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
    # ckpt_dir = "scratch"
    if not os.path.exists(outdir):
        raise Exception("could not find {}".format(outdir))
        
    with open(os.path.join(outdir, "initial_vm_state.json"), "w") as f:
        f.write(json.dumps(vm_initial_state.to_json(), indent=4))
    sim_model.load_from_json_state(vm_initial_state.to_json())#, verbose=True)
    print("Finished sim_mode.load_from_json_state")
    
    checkpoint_fps = []
    
    
    for fn in os.listdir(outdir):
        if (fn.startswith("res") or fn.startswith("vmst_")) and fn.endswith(".json"):
            os.remove(os.path.join(outdir, fn))
    
    ckpt_strnum_nchars = len(str(args.n_checkpoints - 1))
    
    tot_time_elaspsed = 0
    for i in range(args.n_checkpoints):
        ckpt_filename = "res{}.json".format(str(i).zfill(ckpt_strnum_nchars))
        ckpt_fp = os.path.join(outdir, ckpt_filename)
        vmstate_json = sim_model.vm_state_json()
        
        with open(ckpt_fp, "w") as f:
            f.write(json.dumps(vmstate_json, indent=4))
        
        checkpoint_fps.append(ckpt_fp)
        sim_model.run_with_adaptive_tstep(args.ckpt_period, do_time_force_computation=True, verbose=False)
        tot_time_elaspsed += args.ckpt_period
        print("Time elapsed: {}".format(tot_time_elaspsed))
        if sim_model._check_topology_changed():
            print("TOPOLOGY CHANGED")

@dataclass
class ExperimentConfig:
    box_size_rel_x: float
    box_size_rel_y: float
    field_strength: float
    param_A0: float
    param_P0: float
    param_gamma: float
    param_kappa: float
    field_size_multiplier: float
    vertex_friction_gamma: float
    t1_min_edge_len_rel: float
    t1_new_edge_len_rel: float
    init_integrator_dt: float
    dormand_prince_max_error: float
    ckpt_period: int
    n_checkpoints: int



def generate_manifest_experiment_structure_from_battery_config(jobj, root_outdir):
    # root_outdir 
    if not os.path.isdir(root_outdir):
        raise ValueError("root outir '{}' is not a directory".format(root_outdir))
    
    if "vary_shape_param" in jobj["settings"] and jobj["settings"]["vary_shape_param"]:
        shapeparam_vals = jobj["settings"]["shape_param_values"]
        
        assert "param_P0" not in jobj["defaults"]
        assert "param_A0" not in jobj["defaults"]
        
        param_A0 = None
        param_P0 = None
        
        if not isinstance(shapeparam_vals, list):
            raise ValueError("Expected jobj/settings/shape_param_values to be a list")
        # for sp in shapeparam_vals:
    else:
        assert "shape_param_values" not in jobj["settings"]
        shapeparam_vals = None
    
    param_gamma = jobj["defaults"]["param_gamma"]
    param_kappa = jobj["defaults"]["param_kappa"]
    
    
    
    ### Set up range of params
    vals_for_A0_P0 = []
    if shapeparam_vals is not None:
        # Generate a range of A0, P0 values that:
        #  1. have the correct shape parameters specified
        #  2. will have a rest area of 1.0
        for shapep in shapeparam_vals:
            A0_unscaled = 1
            P0_unscaled = shapep * math.sqrt(A0_unscaled)
            
            unscaled_rest_side_length = RegHexagonalModel().find_rest_size_of_hexagon(
                A_0=A0_unscaled,
                P_0=P0_unscaled,
                K=param_kappa,
                Gamma=param_gamma,
            )["rest_side_length"]
            
            unscaled_rest_area = ( (3 * (3**0.5))/2 ) * (unscaled_rest_side_length**2)
            
            # We scale these, so that the cells will have a rest area of 1. this makes things easier to compare.
            A0_scaled = A0_unscaled / unscaled_rest_area
            P0_scaled = P0_unscaled / math.sqrt(unscaled_rest_area)
            
            vals_for_A0_P0.append({
                "A0": A0_scaled,
                "P0": P0_scaled,
            })
    else:
        # Just use what was provided
        vals_for_A0_P0.append({
            "A0": jobj["defaults"]["param_A0"],
            "P0": jobj["defaults"]["param_P0"],
        })
    
    all_experiment_configs = []
    for pair_A0P0 in vals_for_A0_P0:
        # print("Build config")
        all_experiment_configs.append({
            "params": {
                "box_size_rel_x":               jobj["defaults"]["box_size_rel_x"],
                "box_size_rel_y":               jobj["defaults"]["box_size_rel_y"],
                "field_strength":           jobj["defaults"]["field_strength"],
                "param_A0":                 pair_A0P0["A0"],
                "param_P0":                 pair_A0P0["P0"],
                "param_gamma":              jobj["defaults"]["param_gamma"],
                "param_kappa":              jobj["defaults"]["param_kappa"],
                "field_size_multiplier":    jobj["defaults"]["field_size_multiplier"],
                "ckpt_period":              jobj["defaults"]["ckpt_period"],
                "vertex_friction_gamma":    jobj["defaults"]["vertex_friction_gamma"],
                "t1_min_edge_len_rel":      jobj["defaults"]["t1_min_edge_len_rel"],
                "t1_new_edge_len_rel":      jobj["defaults"]["t1_new_edge_len_rel"],
                "init_integrator_dt":       jobj["defaults"]["init_integrator_dt"],
                "dormand_prince_max_error": jobj["defaults"]["dormand_prince_max_error"],
                "n_checkpoints":            jobj["defaults"]["n_checkpoints"],
            },
        })
    
    # Generate a main json file to point to all of the child jsons
    manifest_jobj = {
        "sub_experiments": []
    }
    exp_numstrlen = len(str( len(all_experiment_configs) - 1 ))
    for exp_idx, exp_conf in enumerate(all_experiment_configs):
        exp_out_dirname = "exp_{}".format(str(exp_idx).zfill(exp_numstrlen))
        checkpoints_dirname = "ckpts"
        experiment_conf_fn = "exp_conf.json"
        
        # Make directories
        exp_base_dirpath = os.path.join(root_outdir,exp_out_dirname)
        experiment_output_dirpath = os.path.join(exp_base_dirpath, checkpoints_dirname)
        os.mkdir(exp_base_dirpath)
        os.mkdir(experiment_output_dirpath)
        
        exp_conf_abspath = os.path.join(
            root_outdir,
            exp_out_dirname,
            experiment_conf_fn
        )
        # Write the config file
        with open(exp_conf_abspath, "w") as f:
            print(exp_conf)
            print("Writing experiment config to '{}'".format(exp_conf_abspath))
            json.dump(exp_conf, f, indent=4)
        
        manifest_jobj["sub_experiments"].append({
            "experiment_dirname": exp_out_dirname,
            "conf_path": os.path.join(exp_out_dirname, experiment_conf_fn),
            "ckpts_dir": os.path.join(exp_out_dirname, checkpoints_dirname),
        })
    manifest_abspath = os.path.join(
        root_outdir,
        "manifest_experiments.json",
    )
    with open(manifest_abspath, "w") as f:
        print("\nWriting manifest:  {}".format(manifest_abspath))
        json.dump(manifest_jobj, f, indent=4)
    
    return manifest_jobj

def generate_battery_setup(args):
    if not os.path.exists(args.config_json):
        raise ValueError("'{}' does not exist".format(args.config_json))
    if not os.path.isdir(args.output_dir):
        raise ValueError("'{}' is not a directory".format(args.output_dir))
    
    with open(args.config_json, "r") as f:
        battery_config = generate_manifest_experiment_structure_from_battery_config(
            json.load(f),
            root_outdir=args.output_dir,
        )

def spawn_child(run_dirpath, run_config_json_fp):
    print("\nSpawning child...")
    # print(run_dirpath)
    print("Reading from: ", run_config_json_fp)
    with open(run_config_json_fp, "r") as f:
        conf_jobj = json.load(f)
    
    experiment_config = ExperimentConfig(**conf_jobj["params"])
    print(experiment_config)
    
    run_hexagonal(experiment_config, outdir=run_dirpath)


def run_battery(args):
    if not os.path.exists(args.manifest_json):
        raise ValueError("Could not find '{}'".format(args.manifest_json))
    
    with open(args.manifest_json, 'r') as f:
        manifest_jobj = json.load(f)
    
    root_batterydir = os.path.dirname(args.manifest_json)
    
    # If we're working in the same directory, it may just be a blank string - so use cwd '.'
    if root_batterydir == "":
        root_batterydir = "."
    
    for sub_exp in manifest_jobj["sub_experiments"]:
        # {'experiment_dirname': 'exp_0', 'conf_path': 'exp_0/exp_conf.json', 'ckpts_dir': 'exp_0/ckpts'}
        sub_exp["experiment_dirname"]
        sub_exp["conf_path"]
        conf_json_path = os.path.join(root_batterydir, sub_exp["conf_path"])
        ckpts_absdir = os.path.join(root_batterydir, sub_exp["ckpts_dir"])
        
        spawn_child(
            run_dirpath=ckpts_absdir,
            run_config_json_fp=conf_json_path,
        )
        # print(conf_json_path)
        # print(ckpts_absdir)
    
        
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    subparsers = parser.add_subparsers(
        help="Which command to run - manual single experiment, or a range of tests",
        dest="command",
        required=True,
    )

    single_experiment_parser = subparsers.add_parser("single")
    single_experiment_parser.add_argument("outdir", help="Output dir for the checkpoint files")
    single_experiment_parser.add_argument("--box_size_rel_x", default=30, type=float, help="X size of tissue, in sqrt(equilibrium area).")
    single_experiment_parser.add_argument("--box_size_rel_y", default=30, type=float, help="Y size of tissue, in sqrt(equilibrium area).")
    single_experiment_parser.add_argument("--field_strength", default=1.2, type=float, help="field strength")
    single_experiment_parser.add_argument("--param_A0", default=15, type=float, help="A0 energy term")
    single_experiment_parser.add_argument("--param_P0", default=5, type=float, help="A0 energy term")
    single_experiment_parser.add_argument("--param_gamma", default=1.0, type=float, help="A0 energy term")
    single_experiment_parser.add_argument("--param_kappa", default=1.0, type=float, help="A0 energy term")
    single_experiment_parser.add_argument("--field_size_multiplier", default=0.4, type=float, help="Size of the field")
    single_experiment_parser.add_argument("--ckpt_period", default=0.300, type=float, help="Size of each integration step")
    single_experiment_parser.add_argument("--vertex_friction_gamma", default=0.1 ,type=float, help="Friction coefficient for vertex motion")
    single_experiment_parser.add_argument("--t1_min_edge_len_rel", default=0.2, type=float, help="Minimum edge length for t1 transition, relative to edge length at rest")
    single_experiment_parser.add_argument("--t1_new_edge_len_rel", default=0.21, type=float, help="New edge length after t1 transition, relative to edge length at rest")
    single_experiment_parser.add_argument("--init_integrator_dt", default=0.005, type=float, help="Initial dt to try for the adaptive integrator")
    single_experiment_parser.add_argument("--dormand_prince_max_error", default= 0.02, type=float,help="Max error (between fourth and fifth order runge kutta dormand prince results) to allow for a calcuation of a new vertex position.")
    single_experiment_parser.add_argument("--n_checkpoints", default=500, type=int, help="Number of checkpoints to run")
    
    
    
    battery_generator_parser = subparsers.add_parser("generate_battery_setup")
    battery_generator_parser.add_argument("--config_json", required=True, help="Json file containing the specifications for the experiments to be run")
    battery_generator_parser.add_argument("--output_dir", required=True, help="Directory to put results of experiments")
    
    battery_runner_parser = subparsers.add_parser("run_battery")
    battery_runner_parser.add_argument("--manifest_json", required=True, help="Json file with the manifest of the child experiments to be run. Generated with the generate_battery_setup command.")
    
    
    
    
    
    args = parser.parse_args()
    
    if args.command == "single":
        configs = run_hexagonal(args, outdir=args.outdir)
    elif args.command == "generate_battery_setup":
        generate_battery_setup(args)
    elif args.command == "run_battery":
        run_battery(args)
    else:
        raise ValueError("unknown command '{}'".format(args.command))

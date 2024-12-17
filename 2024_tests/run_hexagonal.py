import numpy as np
import argparse
import json
import os
from simcode.analytic_tools.find_hexagon_rest_area import HexagonalModel, find_hexagon_rest_area

from VMToolkit.config_builder.open.honeycomb_lattice import HoneycombLattice
# from VMToolkit.VM import Tissue, System, ForceCompute, Integrate, Topology, Dump, Simulation, Vec
from VMToolkit.VMAnalysis.utils.HalfEdge import Mesh

from simcode.tissue_builder.hexagonal import HexagonalCellMesh
from simcode.sim_model.sim_model import SimModel
from simcode.sim_model.vm_state import (
    VMState, SimulationSettings, IntegratorSettings, AllForces, CellGroupForces,
    VertexGroupForces, CellAreaForce, CellPerimeterForce, ConstantVertexForce, 
    # VertexGroup,
)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--fpull", default=-0.2,type=float, help="Force to squeeze/stretch")
    parser.add_argument("--box_lx", default=1.0,type=float, help="X size of tissue")
    parser.add_argument("--box_ly", default=1.0,type=float, help="Y size of tissue")
    
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
    
    
    cm = HexagonalCellMesh(
        side_length=rest_side_length*1.1,
        box_lx=args.box_lx,
        box_ly=args.box_ly,
        verbose=True,
    )
    tiss_topology, tiss_init_state = cm.build_vm_state()#verbose=True)
    
    vm_initial_state = VMState(
        tiss_topology=tiss_topology,
        current_state=tiss_init_state,
        sim_settings=SimulationSettings(
            integrator_settings=IntegratorSettings(
                vertex_friction_gamma=0.1,
                step_dt=0.1,
            ),
            force_settings=AllForces(
                cell_forces=CellGroupForces(group_forces={
                    "all": [
                        CellAreaForce(A0=A0_model, kappa=kappa),
                        CellPerimeterForce(gamma=gamma, lam=P0_model*gamma),
                    ]
                }),
                vertex_forces=VertexGroupForces(vertex_group_forces={
                    "all": [
                        ConstantVertexForce(f_x=0.1,f_y=0)
                    ]
            # const_vtx_prop_force_configs.append(const_force_conf_dict)
                    
                }),
            )
        )
    )
    
    
    
    print("MODEL PARAMS")
    print("A0={A0}  P0={P0}  kappa={kappa} gamma={gamma}".format(A0=A0_model,P0=P0_model,kappa=kappa,gamma=gamma))
    print(analytical_predictions)
    print("theoretical rest width={}, height={}".format(theoretical_rest_width, theoretical_rest_height))
    # ##### Running sim
    # tissue  = Tissue()                                               # initialise mesh
    # sim_sys = System(tissue)                                         # base object for the system
    # forces = Force(sim_sys)                                          # handles all types of forces
    # integrators = Integrate(sim_sys, forces, 0)              # handles all integrators
    # topology = Topology(sim_sys, forces)                             # handles all topology changes (T1, division, ingression)
    # dumps = Dump(sim_sys, forces)                                    # handles all data output 
    # simulation = Simulation(sim_sys, integrators, forces, topology)  # simulation object



    # # #################################################################
    # #
    # # Create the initial configuration and read it
    # #
    # # #################################################################
    sim_model = SimModel(verbose=False)
    
    sim_model.load_from_json_state(vm_initial_state.to_json())#, verbose=True)
    # exit()s
    
    # sim_model.configure_forces(P0_model, gamma, kappa, verbose=True)
    
    # sim_model.load_json_obj(
    #     cm.build_vm_state(verbose=True),
    # )
    
    
    # with open("scratch/forgodot.json", "w") as f:
    #     json.dump(sim_model.get_json_state(), f)
    
    # sim_model.configure_integrators(verbose=True)

    step_size =100      # Step counter in terms of time units
    
    ext_forcing_on = []
    checkpoint_fps = []
    # Pulling on the passive system
    
    ckpt_dir = "scratch"
    if not os.path.exists(ckpt_dir):
        raise Exception("could not find {}".format(ckpt_dir))
    
    for fn in os.listdir(ckpt_dir):
        if fn.startswith("res") and fn.endswith(".json"):
            os.remove(os.path.join(ckpt_dir, fn))
    
    N_checkpoints =20
    for i in range(N_checkpoints):
        ckpt_fp = "scratch/res{}.json".format(str(i).zfill(3))
        json_str = sim_model.dump_cpp_json()
        with open(ckpt_fp, "w") as f:
            f.write(json_str)
        
        checkpoint_fps.append(ckpt_fp)
        
        # Will be reocrded next iteration

        if i/N_checkpoints>=0.25:
            pass
            print("Using external forces.. ")
            # sim_model.set_some_forces()
            
            ext_forcing_on.append(True)
        else:
            ext_forcing_on.append(False)
        
        sim_model.run_steps(step_size)#, verbose=True)
        
        # exit()
    
    print("Running analysis on simulation results")
    for ckpt_idx, ckpt_fp in enumerate(checkpoint_fps):

        print(ckpt_fp,end="")

        mesh = Mesh()
        mesh.read(ckpt_fp)

        cell_widths = []
        cell_heights = []
        passive_real_cells= []
        for f in mesh.faces:
            if f.outer == False:
                
                cell_vertices = np.array([v.to_list() for v in f.vertex_coords()],dtype=float)
                
                
                cell_w = cell_vertices[:,0].max() - cell_vertices[:,0].min()
                cell_h = cell_vertices[:,1].max() - cell_vertices[:,1].min()
                cell_widths.append(cell_w)
                cell_heights.append(cell_h)
                
                passive_real_cells.append(f)
        
        # print(mesh.faces[len(mesh.faces)//3].vertex_coords())
        areas = [face.area() for face in passive_real_cells]
        areas = np.array(areas)

        W_mean = np.array(cell_widths).mean()
        H_mean = np.array(cell_heights).mean()

        print("  A_min={:.4f}, A_max={:.4f}, A_mean={:.4f}, W_mean={W_mean:.4f}, H_mean={H_mean:.4f}".format(
            areas.min(),
            areas.max(),
            areas.mean(),
            W_mean=W_mean,
            H_mean=H_mean,
            ))
        if ext_forcing_on[ckpt_idx]:
            strain_x = (W_mean - theoretical_rest_width)/theoretical_rest_width

            strain_y = (H_mean - theoretical_rest_height)/theoretical_rest_height
            print("    strain_x={} strain_y={}, poisson_r={}".format(strain_x, strain_y, -strain_y/strain_x))
        else:
            pass
        
        if ckpt_idx == len(checkpoint_fps) - 1:
            print([v.to_list() for v in passive_real_cells[0].vertex_coords()])
    # print()
    


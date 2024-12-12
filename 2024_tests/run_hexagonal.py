import numpy as np
import argparse
import json
import os
from analytic_tools.find_hexagon_rest_area import HexagonalModel, find_hexagon_rest_area

from VMToolkit.config_builder.open.honeycomb_lattice import HoneycombLattice
from VMToolkit.VM import Tissue, System, Force, Integrate, Topology, Dump, Simulation, Vec
from VMToolkit.VMAnalysis.utils.HalfEdge import Mesh

from tissue_builder.hexagonal import HexagonalCellMesh
from sim_model.sim_model import SimModel


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--fpull", default=-0.2,type=float, help="Force to squeeze/stretch")
    parser.add_argument("--box_lx", default=8.0,type=float, help="X size of tissue")
    parser.add_argument("--box_ly", default=8.0,type=float, help="Y size of tissue")
    
    args = parser.parse_args()
    
    p0_shapefac = 3.6
    A0_model = 1
    P0_model = 0.2
    kappa = 1.0     # area stiffness
    gamma = 0.15    # perimeter stiffness
    
    analytical_predictions = HexagonalModel().find_rest_size_of_hexagon(
        A_0=A0_model,
        P_0=P0_model,
        K=kappa,
        Gamma=gamma,
    )
    
    rest_side_length = analytical_predictions["rest_side_length"]
    theoretical_rest_width = rest_side_length*2
    theoretical_rest_height = rest_side_length*2*np.sqrt(3)/2
    
    if os.path.exists("scratch/example.json"):
        os.remove("scratch/example.json")
    
    cm = HexagonalCellMesh(
        side_length=rest_side_length*0.8,
        box_lx=args.box_lx,
        box_ly=args.box_ly,
    )
    # cm.log_
    cm.set_all_A0(A0_model)
    cm.set_all_P0(P0_model)
    
    with open("scratch/example.json", "w") as f:
        json.dump(cm.build_vm_mesh_obj(verbose=True), f)
    
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
    sim_model = SimModel(verbose=True)
    
    sim_model.configure_forces(P0_model, gamma, kappa, verbose=True)

    sim_model.load_json("scratch/example.json",verbose=True)           # read input configuration
    with open("scratch/forgodot.json", "w") as f:
        json.dump(sim_model.get_json_state(), f)
    
    sim_model.configure_integrators(verbose=True)

    step_size = 500      # Step counter in terms of time units
    
    ext_forcing_on = []
    checkpoint_fps = []
    # Pulling on the passive system
    # for i in range(args.nrun):
    N_checkpoints = 40
    for i in range(N_checkpoints):
        ckpt_fp = "scratch/res{}.json".format(str(i).zfill(3))
        json_str = sim_model.dump_cpp_json()
        with open(ckpt_fp, "w") as f:
            f.write(json_str)
        
        checkpoint_fps.append(ckpt_fp)
        
        # Will be reocrded next iteration

        if i/N_checkpoints>=0.25:
        #if False:
            pass
            print("Using external forces.. ")
            sim_model.set_some_forces()
            # fpull=args.fpull
            # integrators.set_external_force('brownian', 'right', Vec(fpull,0.0))  # pulling on the right-most column of vertices
            # integrators.set_external_force('brownian', 'left', Vec(-fpull,0.0))  # pulling on the left-most column of vertices
            
            ext_forcing_on.append(True)
        else:
            ext_forcing_on.append(False)
        # simulation.run(step_size, topological_change=False)
        sim_model.run_steps(step_size, verbose=True)
    
    print("Running analysis on simulation results")
    for ckpt_idx, ckpt_fp in enumerate(checkpoint_fps):

        print(ckpt_fp,end="")

        mesh = Mesh()
        mesh.read(ckpt_fp)

        cell_widths = []
        cell_heights = []
        passive_real_cells= []
        for f in mesh.faces:
            if f.type=='passive' and f.outer == False:
                #print(f.vertex_coords())
                cell_vertices = np.array([v.to_list() for v in f.vertex_coords()],dtype=float)
                #print(cell_vertices)
                #break
                cell_w = cell_vertices[:,0].max() - cell_vertices[:,0].min()
                cell_h = cell_vertices[:,1].max() - cell_vertices[:,1].min()
                cell_widths.append(cell_w)
                cell_heights.append(cell_h)
                
                passive_real_cells.append(f)
                #break
        print(mesh.faces[len(mesh.faces)//3].vertex_coords())
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
    


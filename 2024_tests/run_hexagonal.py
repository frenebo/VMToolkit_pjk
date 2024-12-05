import numpy as np
import os
from hexagonal_test.analytic_tools.find_hexagon_rest_area import HexagonalModel, find_hexagon_rest_area

from VMToolkit.config_builder.open.make_honeycomb import create_honeycomb_json
from VMToolkit.config_builder.open.honeycomb_lattice import HoneycombLattice
from VMToolkit.VM import Tissue, System, Force, Integrate, Topology, Dump, Simulation, Vec
from VMToolkit.VMAnalysis.utils.HalfEdge import Mesh

def setup_hexagonal_init_mesh(A0_model,P0_model,init_side_length,json_out_fp):
    box_lx = 25.0
    box_ly = 25.0
    h = HoneycombLattice(
        lx=box_lx,
        ly=box_ly,
        a=init_side_length,
    )
    h.build()
    h.minmax()
    h.set_energy_params(A0=A0_model,P0=P0_model)    
    
    # left_box_corners = [
    #     [-1.0, -1],
    #     [-1.0, box_lx+1],
    #     [init_side_length],
    # ]
    print("C point:")
    
    cell_centers = np.array([cell.rc for cell in h.cells])
    leftmost_center_x = cell_centers[:,0].min()
    rightmost_center_x = cell_centers[:,0].max()
    
    left_box_corners = [
        [-box_lx*0.6, -box_ly*0.6],
        [-box_lx*0.6, box_ly*0.6],
        [
            leftmost_center_x + init_side_length*0.01,
            box_ly*0.6,
        ],
        [
            leftmost_center_x + init_side_length*0.01,
            -box_ly*0.6,
        ],
    ]

    right_box_corners = [
        [box_lx*0.6, -box_ly*0.6],
        [box_lx*0.6, box_ly*0.6],
        [
            rightmost_center_x - init_side_length*0.01,
            box_ly*0.6,
        ],
        [
            rightmost_center_x - init_side_length*0.01,
            -box_ly*0.6,
        ],
    ]

    h.set_vertex_type(left_box_corners, "left")
    h.set_vertex_type(right_box_corners, "right")

    h.set_cell_type(left_box_corners, "leftcell")
    h.set_cell_type(right_box_corners, "rightcell")
    left_faces = [face for face in h.cells if face.type=='left']
    right_faces = [face for face in h.cells if face.type=='right']
    print("set types of {} left, {} right cells".format(len(left_faces), len(right_faces)))
    h.json_out(json_out_fp)

if __name__ == "__main__":
    hex_model = HexagonalModel()
    
    p0_shapefac = 3.6
    A0_model = 4
    P0_model = p0_shapefac*np.sqrt(A0_model)
    kappa = 1.0     # area stiffness
    gamma = 1.0     # perimeter stiffness
    
    res = hex_model.find_rest_size_of_hexagon(
        A_0=A0_model,
        P_0=P0_model,
        K=kappa,
        Gamma=gamma,
    )
    rest_side_length = res["rest_side_length"]
    theoretical_rest_width = rest_side_length*2
    theoretical_rest_height = rest_side_length*2*np.sqrt(3)/2
    print("MODEL PARAMS")
    if os.path.exists("scratch/example.json"):
        os.remove("scratch/example.json")
    setup_hexagonal_init_mesh(
        A0_model,
        P0_model,
        init_side_length=rest_side_length*1.05,
        json_out_fp="scratch/example.json",
    )
    
    print("A0={A0}  P0={P0}  kappa={kappa} gamma={gamma}".format(A0=A0_model,P0=P0_model,kappa=kappa,gamma=gamma))
    print(res)
    print("theoretical rest width={}, height={}".format(theoretical_rest_width, theoretical_rest_height))
    ##### Running sim
    tissue  = Tissue()                                               # initialise mesh
    sim_sys = System(tissue)                                         # base object for the system
    forces = Force(sim_sys)                                          # handles all types of forces
    integrators = Integrate(sim_sys, forces, 0)              # handles all integrators
    topology = Topology(sim_sys, forces)                             # handles all topology changes (T1, division, ingression)
    dumps = Dump(sim_sys, forces)                                    # handles all data output 
    simulation = Simulation(sim_sys, integrators, forces, topology)  # simulation object

    # #################################################################
    #
    # Create the initial configuration and read it
    #
    # #################################################################

    sim_sys.read_input("scratch/example.json")           # read input configuration


    # #################################################################
    #
    # Add forces to the system
    #
    # #################################################################

        

    forces.add('area')         # add area force form term E = 0.5*kappa*(A-A0)^2
    forces.add('perimeter')    # add perimeter force term from E = 0.5*gamma*P^2 + lambda*P (maybe -?)

    # Set parameters for each cell type

    lambda_val = P0_model * gamma # @TODO either the sim uses this, or it uses the P0... add a way to force P0 usage in set_params in cpp file?

    for c_type in ['leftcell', 'rightcell','passive']:
        forces.set_params('area', c_type, {'kappa' : kappa})
        forces.set_params('perimeter', c_type,  {'gamma': gamma, "lambda": lambda_val})

    
    # #################################################################
    #
    # Set conditions for the T1 transition
    #
    # #################################################################

    topology.set_params({'min_edge_len': 0.05, 'new_edge_len': 0.055}) 


    # #################################################################
    #
    # Add Brownian integrator that will handle mechanical part
    #
    # #################################################################

    integrators.add('brownian')    

    # #################################################################
    #
    # Simulation starts here
    #
    # #################################################################

    dt = 0.05
    friction_gam = 1.0
    integrators.set_dt(dt) # set time step
    integrators.set_params("brownian", {"gamma": friction_gam})

    step_size = 200      # Step counter in terms of time units
    
    ext_forcing_on = []
    checkpoint_fps = []
    # Pulling on the passive system
    # for i in range(args.nrun):
    N_checkpoints =50
    for i in range(N_checkpoints):
        ckpt_fp = "scratch/res{}.json".format(str(i).zfill(3))
        dumps.dump_mesh(ckpt_fp)
        checkpoint_fps.append(ckpt_fp)
        
        # Will be reocrded next iteration

        if i/N_checkpoints>=0.3:
        #if False:
            pass
            print("Using external forces.. ")
            fpull=0.1
            integrators.set_external_force('brownian', 'right', Vec(fpull,0.0))  # pulling on the right-most column of vertices
            integrators.set_external_force('brownian', 'left', Vec(-fpull,0.0))  # pulling on the left-most column of vertices
    
            ext_forcing_on.append(True)
        else:
            ext_forcing_on.append(False)
        simulation.run(step_size)
    
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
            #print("")
    #### Do stretching test

    


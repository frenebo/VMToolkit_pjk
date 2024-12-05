import numpy as np
from hexagonal_test.analytic_tools.find_hexagon_rest_area import HexagonalModel, find_hexagon_rest_area

from VMToolkit.config_builder.open.make_honeycomb import create_honeycomb_json
from VMToolkit.config_builder.open.honeycomb_lattice import HoneycombLattice
from VMToolkit.VM import Tissue, System, Force, Integrate, Topology, Dump, Simulation
from VMToolkit.VMAnalysis.utils.HalfEdge import Mesh

if __name__ == "__main__":
    hex_model = HexagonalModel()
    
    A0_model = 20.0
    P0_model = 1.0
    kappa = 1.0     # area stiffness
    gamma = 1.0     # perimeter stiffness
    print("MODEL PARAMS")
    print("A0={A0}  P0={P0}  kappa={kappa} gamma={gamma}".format(A0=A0_model,P0=P0_model,kappa=kappa,gamma=gamma))
    res = hex_model.find_rest_size_of_hexagon(
        A_0=A0_model,
        P_0=P0_model,
        K=kappa,
        Gamma=gamma,
    )
    rest_side_length = res["rest_side_length"]
    print(res)
    
    h = HoneycombLattice(
        lx=25.0,
        ly=25.0,
        a=rest_side_length*2.0,
    )
    h.build()
    h.minmax()
    h.set_energy_params(A0=A0_model,P0=P0_model)    
    h.json_out("scratch/example.json")
    
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

    forces.set_params('area', 'passive', {'kappa' : kappa})
    forces.set_params('perimeter', 'passive',  {'gamma': gamma, "lambda": lambda_val})

    
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

    dt = 0.005
    friction_gam = 2.0
    integrators.set_dt(dt) # set time step
    integrators.set_params({"gamma": friction_gam})

    step_size = 10      # Step counter in terms of time units
    
    checkpoint_fps = []
    # Pulling on the passive system
    # for i in range(args.nrun):
    for i in range(10):
        ckpt_fp = "scratch/res{}.json".format(str(i).zfill(3))
        dumps.dump_mesh(ckpt_fp)
        checkpoint_fps.append(ckpt_fp)
        
        # Will be reocrded next iteration
        simulation.run(step_size)

    print("Running analysis on simulation results")
    for ckpt_fp in checkpoint_fps:
        print(ckpt_fp,end="")

        mesh = Mesh()
        mesh.read(ckpt_fp)
        areas = [face.area() for face in mesh.faces]
        areas = np.array(areas)

        print("  Min,max,mean area: {},{},{}".format(areas.min(), areas.max(),areas.mean()))


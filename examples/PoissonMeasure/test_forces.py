import os
import argparse
import math

from VMToolkit.config_builder.open.make_honeycomb import create_honeycomb_json

from VMToolkit import VM as vm

def test_forces(scratch_dirpath):
    os.makedirs(scratch_dirpath)
    side_len_rest = 1.0
    
    ## Setup initial condition file
    honeycomb_pth = os.path.join(scratch_dirpath, "start_honeycomb.json")
    create_honeycomb_json(honeycomb_pth, honeycomb_a_len=side_len_rest)
    
    
    tissue  = vm.Tissue()                                               # initialise mesh
    sim_sys = vm.System(tissue)                                         # base object for the system
    forces = vm.Force(sim_sys)                                          # handles all types of forces
    integrators = vm.Integrate(sim_sys, forces, seed)                   # handles all integrators
    topology = vm.Topology(sim_sys, forces)                             # handles all topology changes (T1, division, ingression)
    dumps = vm.Dump(sim_sys, forces)                                    # handles all data output 
    simulation = vm.Simulation(sim_sys, integrators, forces, topology)  # simulation object


    sim_sys.read_input(honeycomb_pth)           # read input configuration

    forces.add('area')         # add area force form term E = 0.5*kappa*(A-A0)^2
    forces.add('perimeter')    # add perimeter force term from E = 0.5*gamma*P^2 + lambda*P (maybe -?)
    ###### How do we set gamma and lambda?
    # Making it fit the (P-P0)^2 form, we complete the square to get:
    # 0.5*gamma*P^2+lambda*P = 0.5*g*(P+l/g)^2 + some constant.....
    # So, g = k_perimeter
    # l/g = P_0 -> l = g * P_0
    # And then 
    k_perimeter = 1.0
    perimeter_rest = side_len_rest * 6
    gamma = k_perimeter
    lam = gamma * perimeter_rest
    
    
    area_rest = side_len_rest ** 2 + 3*math.sqrt(3)/2

    # Set parameters for each cell type
    forces.set_params('area', 'passive' , {'kappa' : kappa})
    forces.set_params('perimeter', 'passive' ,  {'gamma': gamma, 'lambda': lam})
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    
    parser.add_argument("--scratch_dir", default="scratch", help="Name of scratch directory to make and create temporary files")
    args = parser.parse_args()
    
    if os.path.exists(args.scratch_dir) and not os.path.isdir(args.scratch_dir):
        raise ValueError("Cannot use {} as scratch dir, it exists but is not a directory".format(args.scratch_dir))
    
    test_forces(scratch_dirpath=args.scratch_dir)
    
    
    
    

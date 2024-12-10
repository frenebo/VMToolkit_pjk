import json
from VMToolkit.VM import Tissue, System, Force, Integrate, Topology, Dump, Simulation, Vec

class VMWrapper:
    def __init__(self):
        raise NotImplementedError()

class SimModel:
    def __init__(self, verbose):
        self._sim_sys = None
        self._forces = None
        self._integrators = None
        self._dumps = None
        self._simulation = None
        
        
        self._initialize_cpp(verbose=verbose)
        
        self.forces_configured = False
        self._tissue_loaded = False
        pass
    
    def configure_integrators(self,verbose=False):
        if not self._tissue_loaded:
            raise Exception("Cannot configure integrator before tissue has been loaded" +
                " - currently c++ code relies on tissue types being there before brownian integrator")
        if verbose:
            print("Adding brownian integrator")
        self._integrators.add('brownian')    


        dt = 0.08
        friction_gam = 1.0
        if verbose:
            print("Setting dt={}, friction_gamma={}".format(dt, friction_gam))
        self._integrators.set_dt(dt) # set time step
        self._integrators.set_params("brownian", {"gamma": friction_gam})
        if verbose:
            print("Done configuring integrators")
    
    def dump_json(self, json_out_fp):
        self._dumps.dump_mesh(json_out_fp)

    def configure_forces(self, P0_model, gamma, kappa, verbose=False):
        # if not self._cell_types_configured:
        #     raise Exception("Cell type not configured - forces won't work without cell types preestablished")
        # if not self._cell_types_
        if verbose:
            print("CONFIGURING FORCES==================")
            self._sim_sys.log_debug_stats()
        self._forces.add('area')         # add area force form term E = 0.5*kappa*(A-A0)^2
        self._forces.add('perimeter')    # add perimeter force term from E = 0.5*gamma*P^2 + lambda*P (maybe -?)
        if verbose:
            print("Added area and perimeter forces")
            self._sim_sys.log_debug_stats()
        # Set parameters for each cell type

        lambda_val = P0_model * gamma # @TODO either the sim uses this, or it uses the P0... add a way to force P0 usage in set_params in cpp file?

        for c_type in ['passive']:
            self._sim_sys.add_cell_type(c_type)
            self._forces.set_params('area', c_type, {'kappa' : kappa})
            self._forces.set_params('perimeter', c_type,  {'gamma': gamma, "lambda": lambda_val})
        
        self.forces_configured = True
    
    def load_json(self, json_fp, verbose=False):
        if not self.forces_configured:
            raise Exception("need to configure forces before loading json")
        
        if verbose:
            print("Loading JSON")
        self._sim_sys.read_input(json_fp)           # read input configuration
        if verbose:
            print("Done loading JSON")
        
        self._tissue_loaded = True

    
    def _initialize_cpp(self, verbose=False):
        ##### Running sim
        tissue  = Tissue()                                               # initialise mesh
        sim_sys = System(tissue)                                         # base object for the system
        if verbose:
            print("constructing force")
            sim_sys.log_debug_stats()
        forces = Force(sim_sys)                                          # handles all types of forces
        if verbose:
            print("constructing integrators")
            sim_sys.log_debug_stats()
        integrators = Integrate(sim_sys, forces, 0)              # handles all integrators
        if verbose:
            print("constructing topology")
            sim_sys.log_debug_stats()
        topology = Topology(sim_sys, forces)                             # handles all topology changes (T1, division, ingression)
        if verbose:
            print("constructing dumps")
            sim_sys.log_debug_stats()
        dumps = Dump(sim_sys, forces)                                    # handles all data output 
        if verbose:
            print("constructing simulation")
            sim_sys.log_debug_stats()
        simulation = Simulation(sim_sys, integrators, forces, topology)  # simulation object                      # handles all data output 
        
        
        """
        IMPORTANT! All these need to be held as CLASS variables python, otherwise data will be lost and bugs will appear.
        """
        self._tissue = tissue # 
        self._sim_sys = sim_sys
        self._forces = forces
        self._integrators = integrators
        self._topology = topology
        self._dumps = dumps
        self._simulation = simulation
        
        if verbose:
            print("Finished intializing cpp...")
            sim_sys.log_debug_stats()
        
    def load_from_hexmesh(self):
        # @TODO find out whether the order of loading force, etc., matters
        raise NotImplementedError()
    
    def run_steps(self, n_steps, verbose=False):
        if not self.forces_configured:
            raise Exception("Need to configure force before running")
        
        if verbose:
            print("About to run simulation for {} steps".format(n_steps))
        
        self._simulation.run(n_steps, topological_change=False,)
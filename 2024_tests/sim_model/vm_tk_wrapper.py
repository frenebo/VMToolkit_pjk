import json
from VMToolkit.VM import Tissue, System, ForceCompute, Integrate, Topology, Dump, Simulation, Vec

# class 
# class CppJson
class CppJsonTissueBuilder:
    @staticmethod
    def build_json_obj_for_vmtoolkit_loader(vm_state):
        
        vertex_id_to_cpp_idx = {}
        cell_id_to_cpp_face_idx = {}
        
        
        # Make 'blank' lists with the length we expect (length = max index plus one)
        cpp_vertex_objs = [None] * max(list(vertex_id_to_cpp_idx.values()) + 1)
        cpp_face_objs = [None] * max(list(face_id_to_cpp_idx.values()) + 1)
        # for 
        for cell_id, cell_top in vm_state.topology().cells().items():
            face_index = cell_id_to_cpp_face_idx[cell_id]
            
            cpp_face_objs[face_index] = {
                "id": face_index,
                "vertices": cell_top.vertex_ids(),
                "nsides": len(cell_top.vertex_ids())
                "outer": is_outer_face,
            }
        
        for vertex_id, vertex_top in vm_state.topology().vertices().items():
            vertex_index = vertex_id_to_cpp_idx[vertex_id]
            
            v_geometry = vm_state.current_state().geometry().vertices()[vertex_id]
            vtx_x = v_geometry.x()
            vtx_y = v_geometry.y()
            
            cpp_vertex_objs[vertex_index] = {
                "id": vertex_index, # This is just so the vertex in the model knows its own index
                "r": [vtx_x, vtx_y],
                "boundary": cell_top.is_boundary(),
                "coorination":num_faces,
                "erased": False,
                "neighbours": vtx_neighbours,
            }
            
        for i, el in enumerate(cpp_vertex_objs):
            if el is None:
                raise Exception("No vertex corresponding to index {}".format(i))
            
        for i, el in enumerate(cpp_face_objs):
            if el is None:
                raise Exception("No face corresponding to index {}".format(i))
        
        return {
            "mesh": {
                # "time_step": 0,
                "vertices": cpp_vertex_objs,
                "faces": cpp_face_objs,
            }
        }

class VMToolkitWrapper:
    def __init__(self, verbose=False):
        self._sim_sys = None
        self._forces = None
        self._integrators = None
        self._dumps = None
        self._simulation = None
        
        # self._state_sinitialized = False
        self._tissue_initialized = False
        
        
        self._initialize_cpp(verbose=verbose)
        
        self.forces_configured = False
        self._tissue_loaded = False
    
    def initialize_from_vm_state(self, vm_state):
        if self._tissue_initialized:
            raise Exception("Cannnot initialize from vm state, when tissue has already been initialized!")
            
        
        cpp_state_json_str = json.dumps(CppJsonTissueBuilder.build_json_obj_for_vmtoolkit_loader(vm_state))
        self._sims_sys.read_input_from_jsonstring(cpp_state_json_str)
        
        self._tissue_initialized = True
        # Build JSON stuff to load in
        # Build the geometry
        
        # Configure forces
        # Configure integrators
        # Comfigure tisues

    
    def _initialize_cpp(self, verbose=False):
        ##### Running sim
        tissue  = Tissue()                                               # initialise mesh
        sim_sys = System(tissue)                                         # base object for the system
        
        if verbose:
            print("constructing force")
            sim_sys.log_debug_stats()
        forces = ForceCompute(sim_sys)  
                                                # handles all types of forces
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
        IMPORTANT! All these need to be held as CLASS variables python, otherwise their memory will be freed,
        and they'll try to talk to each other, and just read whatever is in that memory now. Will cause problems
        """
        self._tissue = tissue
        self._sim_sys = sim_sys
        self._forces = forces
        self._integrators = integrators
        self._topology = topology
        self._dumps = dumps
        self._simulation = simulation
        
        if verbose:
            print("Finished intializing cpp...")
            sim_sys.log_debug_stats()
        
    def dump_cpp_json(self):
        return self._dumps.mesh_to_jsonstr()

    def configure_forces(
        self,
        P0_model,
        gamma,
        kappa,
        verbose=False,
        ):
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

        # lambda_val = P0_model * gamma # @TODO either the sim uses this, or it uses the P0... add a way to force P0 usage in set_params in cpp file?
        self._default_gamma = gamma
        self._default_lambda_val = P0_model * gamma
        self._default_kappa = kappa


        self.forces_configured = True

    def configure_integrators(self, dt, friction_gam, verbose=False):
        # if not self._tissue_loaded:
        #     raise Exception("Cannot configure integrator before tissue has been loaded" +
        #         " - currently c++ code relies on tissue types being there before brownian integrator")
        if verbose:
            print("Adding brownian integrator")
        self._integrators.add('brownian')    


        dt = 0.3
        friction_gam = 1.0
        if verbose:
            print("Setting dt={}, friction_gamma={}".format(dt, friction_gam))
        self._integrators.set_dt(dt) # set time step
        self._integrators.set_params("brownian", {"gamma": friction_gam})
        if verbose:
            print("Done configuring integrators")
    
    def run_steps(
        self,
        n_steps,
        verbose=False,
    ):
        if not self.forces_configured:
            raise Exception("Need to configure force before running")
        
        if verbose:
            print("About to run simulation for {} steps".format(n_steps))
        
        self._simulation.run(n_steps, topological_change=False,)
        
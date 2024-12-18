import json
from VMToolkit.VM import Tissue
from VMToolkit.VM import System
from VMToolkit.VM import ForceCompute, Integrate
from VMToolkit.VM import Topology, Dump, Simulation, Vec

from ..vm_state import CellAreaForce, CellPerimeterForce, ConstantVertexForce
# from .
# class 
# class CppJson
class CppJsonTissueBuilder:
    @staticmethod
    def build_json_obj_for_vmtoolkit_loader(vm_state, vertex_id_to_cpp_idx, cell_id_to_cpp_face_idx):
        
        # raise NotImplementedError("need to convert bween vvertex id and cpp indices")
        
        
        # Make 'blank' lists with the length we expect (length = max index plus one)
        cpp_vertex_objs = [None] * (1+max(list(vertex_id_to_cpp_idx.values())))
        cpp_face_objs = [None] * (1+max(list(cell_id_to_cpp_face_idx.values())))
        
        # vtx_num_faces = {}
        # vtx_neighbours = {}
        # for 
        for cell_id, cell_top in vm_state.topology().cells().items():
            face_index = cell_id_to_cpp_face_idx[cell_id]
            
            cpp_vids = [ vertex_id_to_cpp_idx[vid] for vid in cell_top.vertex_ids() ]
            
            # for vid in cell_top.vertex_ids():
            #     if vid not in vtx_num_faces:
            #         vtx_num_faces[vid] = 0
                
            #     vtx_num_faces[vid] += 1
            
            cpp_face_objs[face_index] = {
                "id": face_index,
                "vertices": cpp_vids,
                "nsides": len(cell_top.vertex_ids()),
                "outer": cell_top.is_outer(),
            }
        
        for vertex_id, vertex_top in vm_state.topology().vertices().items():
            vertex_index = vertex_id_to_cpp_idx[vertex_id]
            
            v_geometry = vm_state.current_state().geometry().vertices()[vertex_id]
            vtx_x = v_geometry.x()
            vtx_y = v_geometry.y()
            
            # num_faces = 
            
            cpp_vertex_objs[vertex_index] = {
                "id": vertex_index, # This is just so the vertex in the model knows its own index
                "r": [vtx_x, vtx_y],
                "boundary": vertex_top.is_boundary(),
                # "coorination": vtx_num_faces[vertex_id],
                "erased": False,
                "constraint": "none",
                # "neighbours": vtx_neighbours,
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
    @staticmethod
    def _make_ids_to_cpp_index_maps(vm_state):
        cell_ids_to_idx = {}
        for cell_idx, cell_id in enumerate(vm_state.topology().cells().keys()):
            cell_ids_to_idx[cell_id] = cell_idx
        
        vtx_ids_to_idx = {}
        for v_idx, vertex_id in enumerate(vm_state.topology().vertices().keys()):
            vtx_ids_to_idx[vertex_id] = v_idx
        
        return {
            "cell_ids_to_idx": cell_ids_to_idx,
            "vtx_ids_to_idx": vtx_ids_to_idx,
            # "vtx_indic"
        }
    def __init__(self, verbose=False):
        self._sim_sys = None
        self._forces = None
        self._integrators = None
        self._dumps = None
        self._simulation = None
        
        # self._state_sinitialized = False
        self._tissue_initialized = False
        self._ids_to_cpp_index_maps = None
        
        
        self._initialize_cpp(verbose=verbose)
        
        self.forces_configured = False
        self._tissue_loaded = False
    
    def initialize_from_vm_state(self, vm_state, verbose=False):
        if self._tissue_initialized:
            raise Exception("Cannnot initialize from vm state, when tissue has already been initialized!")
            
        self._ids_to_cpp_index_maps = VMToolkitWrapper._make_ids_to_cpp_index_maps(vm_state)
        cpp_state_json_str = json.dumps(CppJsonTissueBuilder.build_json_obj_for_vmtoolkit_loader(
            vm_state=vm_state,
            # cellvvver
            vertex_id_to_cpp_idx=self._ids_to_cpp_index_maps["vtx_ids_to_idx"],
            cell_id_to_cpp_face_idx=self._ids_to_cpp_index_maps["cell_ids_to_idx"],
        ))
        # print(cpp_state_json_str)
        self._sim_sys.read_input_from_jsonstring(cpp_state_json_str, verbose=verbose)
        
        self._configure_integrators(
            dt=vm_state.sim_settings().integrator_settings().step_dt(),
            friction_gam=vm_state.sim_settings().integrator_settings().vtx_friction_gamma(),
        )
        self._configure_forces(vm_state, verbose=verbose)
        
        self._tissue_initialized = True
    
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
    
    def _configure_forces(self, vm_state, verbose=False):
        print("Configuring forces")
        cell_groups = vm_state.current_state().cell_groups()
        vertex_groups = vm_state.current_state().vertex_groups()
        
        force_sets = vm_state.sim_settings().force_settings()
        
        
        # cell_area_forces_enabled = []
        area_forces_by_cell = {}
        perimeter_forces_by_cell = {}
        
        
        for cell_grp_name, cell_forces in force_sets.cell_group_forces().force_lists_by_group().items():
            cell_ids_within_group = cell_groups[cell_grp_name].cell_ids()
            for c_force in cell_forces:
                if isinstance(c_force, CellAreaForce):
                    for cell_id in cell_ids_within_group:
                        if cell_id in area_forces_by_cell:
                            raise ValueError("Conflicting cell area forces for cell id {}".format(cell_id))
                        area_forces_by_cell[cell_id] = {
                            "A0": c_force.A0(),
                            "kappa": c_force.kappa(),
                        }
                elif isinstance(c_force, CellPerimeterForce):
                    for cell_id in cell_ids_within_group:
                        if cell_id in perimeter_forces_by_cell:
                            raise  ValueError("Conflicting perimeter forces for cell id {}".format(cell_id))
                        perimeter_forces_by_cell[cell_id] = {
                            "gamma": c_force.gamma(),
                            "lambda": c_force.lam(),
                            # "A0": c_force.A0
                        }
                else:
                    print(c_force)
                    raise Exception("Unkonwo cell force type... {}".format(c_force))
        
        constant_vtx_forces_by_vertex = {}
        for vtx_grp_name, vtx_forces in force_sets.vertex_group_forces().force_lists_by_group().items():
            # print("")
            vtx_ids_within_group = vertex_groups[vtx_grp_name].vertex_ids()
            for v_force in vtx_forces:
                if isinstance(v_force, ConstantVertexForce):
                    for vtx_id in vtx_ids_within_group:
                        if vtx_id in constant_vtx_forces_by_vertex:
                            raise ValueError("Conflicting constant vertex force for cell id {}".format(vtx_id))
                        constant_vtx_forces_by_vertex[vtx_id] = {
                            "f_x": v_force.f_x(),
                            "f_y": v_force.f_y(),
                        }
                else:
                    print(v_force)
                    raise Exception("Unknown vertex force type... {}".format(v_force))
            
        
        ##### Set up these parameters in the model
        self._forces.add('area')         # add area force form term E = 0.5*kappa*(A-A0)^2
        area_force_cell_indices = []
        area_force_configs = []
        for cell_id, area_conf_dict in area_forces_by_cell.items():
            cell_idx = self._ids_to_cpp_index_maps["cell_ids_to_idx"][cell_id]
            
            area_force_cell_indices.append(cell_idx)
            area_force_configs.append(area_conf_dict)
        
        self._forces.set_face_params_facewise(
            "area", area_force_cell_indices, area_force_configs, verbose
        )
        
        
        ### Set up perimeter forces
        self._forces.add('perimeter')    # add perimeter force term from E = 0.5*gamma*P^2 + lambda*P (maybe -?)
                    
                    
        perim_force_cell_indices = []
        perim_force_configs = []
        for cell_id, perim_conf_dict in perimeter_forces_by_cell.items():
            cell_idx = self._ids_to_cpp_index_maps["cell_ids_to_idx"][cell_id]
            perim_force_cell_indices.append(cell_idx)
            perim_force_configs.append(perim_conf_dict)
        
        self._forces.set_face_params_facewise(
            "perimeter", perim_force_cell_indices, perim_force_configs, verbose
        )
        
        self._forces.add("const_vertex_propulsion")
        # constant_vtx_forces_by_vertex = {}
        const_vtx_prop_force_vertex_indices = []
        const_vtx_prop_force_configs = []
        for vtx_id, const_force_conf_dict in constant_vtx_forces_by_vertex.items():
            vtx_idx = self._ids_to_cpp_index_maps["vtx_ids_to_idx"][vtx_id]
            const_vtx_prop_force_vertex_indices.append(vtx_idx)
            const_vtx_prop_force_configs.append(const_force_conf_dict)
        
        self._forces.set_vertex_params_vertexwise(
            "const_vertex_propulsion", const_vtx_prop_force_vertex_indices, const_vtx_prop_force_configs, verbose
        )
            
            
        # self._forces.set_face_params_facewise(
        #     "perimeter", all_face_ids, [{'gamma': self._default_gamma, "lambda": self._default_lambda_val} for i in range(len(all_face_ids))]
        # )
            # print(vtx_grp_name)
            # print(vtx_forces)
            
            # print(vtx_forces)

    # def configure_forces(
    #     self,
    #     P0_model,
    #     gamma,
    #     kappa,
    #     verbose=False,
    #     ):
    #     # if not self._cell_types_configured:
    #     #     raise Exception("Cell type not configured - forces won't work without cell types preestablished")
    #     # if not self._cell_types_
    #     if verbose:
    #         print("CONFIGURING FORCES==================")
    #         self._sim_sys.log_debug_stats()
        
    #     if verbose:
    #         print("Added area and perimeter forces")
    #         self._sim_sys.log_debug_stats()
    #     # Set parameters for each cell type

    #     # lambda_val = P0_model * gamma # @TODO either the sim uses this, or it uses the P0... add a way to force P0 usage in set_params in cpp file?
    #     self._default_gamma = gamma
    #     self._default_lambda_val = P0_model * gamma
    #     self._default_kappa = kappa


    #     self.forces_configured = True

    def _configure_integrators(self, dt, friction_gam, verbose=False):
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
        
        if verbose:
            print("About to run simulation for {} steps".format(n_steps))
        
        self._simulation.run(n_steps, topological_change=False, verbose=verbose)
        
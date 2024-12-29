import json
from VMToolkit.VM import Tissue
from VMToolkit.VM import System
from VMToolkit.VM import ForceCompute, Integrate
from VMToolkit.VM import Topology, Dump, Simulation, Vec

from ..vm_state import (
    CellAreaForce,
    CellPerimeterForce,
    ConstantVertexForce,
    VertexForce,
    CellForce,
    ElectricForceOnCellBoundary,
    EFieldSpecConstantPolygonRegion,
    PolygonSpec,
)
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
        }
    
    def __init__(self, verbose=False):
        self._sim_sys = None
        self._forces = None
        self._integrators = None
        self._dumps = None
        self._simulation = None
        
        self._tissue_initialized = False
        self._ids_to_cpp_index_maps = None
        
        self._last_vm_state = None
        
        
        self._initialize_cpp(verbose=verbose)
    
    def start_force_compute_timers(self, verbose=False):
        self._forces.start_force_compute_timers(verbose=verbose)
        
    def get_force_compute_timers_millis(self,verbose=False):
        return self._forces.get_force_compute_timers_millis(verbose=verbose)
    
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
        self._last_vm_state = vm_state.from_json(vm_state.to_json())
    
    def run_steps(
        self,
        n_steps,
        verbose=False,
    ):
        if not self._tissue_initialized:
            raise Exception("No tissue has been loaded...")
        if verbose:
            print("About to run simulation for {} steps".format(n_steps))
        
        self._simulation.run(n_steps, topological_change=False, verbose=verbose)
        # self._update_vm_state_from_cpp_vm()
        
    def _update_vm_state_from_cpp_vm(self):
        # All that we need to look at (for now) should be the coordinates of vertices - 
        # we are asssuming the topology isn't changing.
        
        # new_vmstate = vm_state.from_json(self._last_vm_state.to_json())
        
        # Get all of the vertex positions
        # new_
        raise NotImplementedError()
        # print("NOT IMP")
    
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
        topology = Topology(sim_sys)                             # handles all topology changes (T1, division, ingression)
        
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
    
    def _add_cell_area_force(self, force_id, force_spec, cell_indices, verbose=False):
        assert isinstance(force_spec, CellAreaForce)
        
        self._forces.add_force(force_id, 'area')         # add area force form term E = 0.5*kappa*(A-A0)^2
        
        area_force_configs = []
        for i in range(len(cell_indices)):
            area_force_configs.append({
                "A0": force_spec.A0(),
                "kappa": force_spec.kappa(),
            })
            
        self._forces.set_face_params_facewise(
            force_id, cell_indices, area_force_configs, verbose
        )
    
    def _add_cell_perimeter_force(self, force_id, force_spec, cell_indices, verbose=False):
        assert isinstance(force_spec, CellPerimeterForce)
        
        self._forces.add_force(force_id, "perimeter")
        
        perimeter_force_configs = []
        for i in range(len(cell_indices)):
            perimeter_force_configs.append({
                "gamma": force_spec.gamma(),
                "lambda": force_spec.lam(),
            })
        
        self._forces.set_face_params_facewise(
            force_id, cell_indices, perimeter_force_configs, verbose
        )
        
    def _add_cell_electric_boundary_force(self, force_id, force_spec, cell_indices, verbose=False):
        assert isinstance(force_spec, ElectricForceOnCellBoundary)
        
        self._forces.add_force(force_id, "force_efield_on_cell_boundary")
        
        field_spec = force_spec.electric_field_spec()
        if not isinstance(field_spec, EFieldSpecConstantPolygonRegion):
            raise ValueError("Unimplemented type of field spec: {}".format(field_spec))
            
        zone_bounds = field_spec.zone_bounds()
        if not isinstance(zone_bounds, PolygonSpec):
            raise ValueError("Unimplemented zone type: {}".format(zone_bounds))
        
        num_params = {
            "E_x": field_spec.E_x(),
            "E_y": field_spec.E_y(),
            "n_polygon_vertices": len(zone_bounds.polygon_vertices()),
        }
        str_params = {
            "field_type": "constant",
            "region_type": "polygon",
        }
        for vert_idx, (vert_x, vert_y) in enumerate(zone_bounds.polygon_vertices()):
            num_params["poly_x" + str(vert_idx)] = vert_x
            num_params["poly_y" + str(vert_idx)] = vert_y
        
        self._forces.set_global_params(force_id, 
            num_params=num_params,
            str_params=str_params,
        )
        
        fids = list(cell_indices)
        face_params = [{} for i in range(len(fids))]
        
        self._forces.set_face_params_facewise(force_id, 
            fids,
            face_params,
            verbose,
        )
    
    def _add_constant_vertex_force(self, force_id, force_spec, vertex_indices, verbose=False):
        assert isinstance(force_spec, ConstantVertexForce)
        
        self._forces.add_force(force_id, "const_vertex_propulsion")
        
        const_vtx_force_configs = []
        for i in range(len(vertex_indices)):
            const_vtx_force_configs.append({
                "f_x": force_spec.f_x(),
                "f_y": force_spec.f_y(),
            })
        
        self._forces.set_vertex_params_vertexwise(
            force_id, vertex_indices, const_vtx_force_configs, verbose
        )
    
    def _add_new_vertex_force(self, force_id, force_spec, vertex_indices_force, verbose=False):
        assert isinstance(force_spec, VertexForce)
        
        if isinstance(force_spec, ConstantVertexForce):
            self._add_constant_vertex_force(force_id, force_spec, vertex_indices_force, verbose=verbose)
        else:
            raise ValueError("Unknown VERTEX force spec type for force id {} ... {}".format(force_id, force_spec))
    
    def _add_new_cell_force(self, force_id, force_spec, cell_indices_force, verbose=False):
        assert isinstance(force_spec, CellForce)
        
        if isinstance(force_spec,  CellAreaForce):
            self._add_cell_area_force(force_id, force_spec, cell_indices_force, verbose=verbose)
        elif isinstance(force_spec, CellPerimeterForce):
            self._add_cell_perimeter_force(force_id, force_spec, cell_indices_force, verbose=verbose)
        elif isinstance(force_spec, ElectricForceOnCellBoundary):
            self._add_cell_electric_boundary_force(force_id, force_spec, cell_indices_force, verbose=verbose)
        else:
            raise ValueError("Unknown CELL force spec type for force id {} ... {}".format(force_id, force_spec))
        
    
    def _configure_forces(self, vm_state, verbose=False):
        if verbose:
            print("Configuring forces")
        
        cell_groups = vm_state.current_state().cell_groups()
        vertex_groups = vm_state.current_state().vertex_groups()
        tiss_forces = vm_state.current_state().forces()
        
        for force_id, force_spec in tiss_forces.items():
            if isinstance(force_spec, VertexForce):
                ##### Find the vertex indicess that this force should be applied to
                vertex_ids_affected = []
                
                for v_grp_id, v_group in vertex_groups.items():
                    if force_id in v_group.force_ids():
                        vertex_ids_affected += [vid for vid in v_group.vertex_ids() if vid not in vertex_ids_affected]
                
                vertex_indices_for_force = [self._ids_to_cpp_index_maps["vtx_ids_to_idx"][vid] for vid in vertex_ids_affected]
                
                # if verboses
                if verbose:
                    print("Configuring force '{}', for vertex indices '{}' - {}".format(force_id, vertex_indices_for_force, force_spec))
                
                self._add_new_vertex_force(force_id, force_spec, vertex_indices_for_force, verbose=verbose)
            elif isinstance(force_spec, CellForce):
                ##### Find the cell indices that this force shoud be applied to
                cell_ids_affected = []
                
                for c_grp_id, c_group in cell_groups.items():
                    if force_id in c_group.force_ids():
                        cell_ids_affected += [cid for cid in c_group.cell_ids() if cid not in cell_ids_affected]
                    
                cell_indices_for_force = [self._ids_to_cpp_index_maps["cell_ids_to_idx"][cid] for cid in cell_ids_affected]
                
                if verbose:
                    print("Configuring force '{}', for vertex indices '{}' - {}".format(force_id, cell_indices_for_force, force_spec))
                
                self._add_new_cell_force(force_id, force_spec, cell_indices_for_force, verbose=verbose)
            else:
                raise ValueError("Expecteed force '{}' to  either be instance of CellForce or VertexForce! {}".format(force_id, force_spec))
        
        

    def _configure_integrators(self, dt, friction_gam, verbose=False):
        if verbose:
            print("Adding runge_kutta integrator")
        self._integrators.add('runge_kutta')
        # if verbose:
        #     print("Adding brownian integrator")
        # self._integrators.add('brownian')


        # dt = 0.3
        # friction_gam = 1.0
        if verbose:
            print("Setting dt={}, friction_gamma={}".format(dt, friction_gam))
        self._integrators.set_dt(dt) # set time step
        self._integrators.set_params("runge_kutta", {"gamma": friction_gam})
        # self._integrators.set_params("brownian", {"gamma": friction_gam})
        if verbose:
            print("Done configuring integrators")
    
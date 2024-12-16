

# class CellTopology:
#     @staticmethod
#     def from_josn(?)

class CellTopology:
    @staticmethod
    def from_json(self, jobj):
        vertex_ids = jobj["vertex_ids"]
        return CellTopology(
            vertices=vertex_ids,
        )
        
    def __init__(self, vertices):
        self._vertex_ids=vertices
    
    def to_json(self):
        return {
            "vertex_ids": self._vertex_ids,
        }

class VertexTopology:
    @staticmethod
    def from_json(self, jobj):
        return VertexTopology()
    
    def __init__(self):
        pass
    
    def to_json(self):
        return {}

class TissueTopology:
    @staticmethod
    def from_json(self, jobj):
        cell_topologies_parsed = {}
        
        for cell_id, cell_top_json in jobj["cells"]:
            cell_topologies_parsed[cell_id] = CellTopology.from_json(cell_top_json)
        
        vertex_topologies_parsed = {}
        for vtx_id, vtx_top_json in jobj["vertices"]:
            vertex_topologies_parsed[vtx_id] = VertexTopology.from_json(vtx_top_json)
        
        # @TODO turn these assertions into checks that give useful exception types and messages that can be handled in a meaningful way
        # assert isinstance(jobj, dict)
        # assert "cells" in jobj
        # assert isinstance(jobj["cells"], dict)
        
        # assert "vertex_ids" in jobj
        # assert isinstance(jobj["vertex_ids"], list)
        
        # for c in jobj["cells"]:
        #     assert isinstance(c, dict)
        #     # assert
        #     assert "verti"
        return TissueTopology(
            cell_topologies=cell_topologies_parsed,
            vertex_topologies=vertex_topologies_parsed
        )
        
    def __init__(self, cell_topologies, vertex_topologies):
        self._cell_topologies = cell_topologies
        self._vertex_topologies = vertex_topologies
    
    
    def to_json(self):
        json_cell_tops = {}
        for cell_id, cell_top in self._cell_topologies:
            json_cell_tops[cell_id] = cell_top.to_json()
        
        json_vtx_tops = {}
        for vtx_id, vtx_top in self._vertex_topologies:
            json_vtx_tops[vtx_id] = vtx_top.to_json()
        return {
            "cells": json_cell_tops,
            "vertex_ids": json_vtx_tops,
        }




class TissueState:
    @staticmethod
    def from_json(jobj):
        return TissueState(
            vertex_geometries=jobj["geometry"]["vertices"],
            vertex_groups=jobj["vertexgroups"],
            cell_groups=jobj["cellgroups"],
        )
        
    def __init__(self, vertex_geometries, vertex_groups, cell_groups):
        self._vertex_geometries = vertex_geometries
        self._vertex_groups = vertex_groups
        self._cell_groups = cell_groups
    
    def _geometry_to_json(self):
        return {
            "vertices": self._vertex_geometries,
        }
    
    def _vertex_groups_to_json(self):
        return self._vertex_groups
    
    def _cell_groups_to_json(self):
        return self._cell_groups
    
    def to_json(self):
        """ Will return something like this:
        {
                "geometry": {
                    "vertices": {
                        "0": {
                            "x": 0.0,
                            "y": 1.0,
                        },
                        etc....
                    }
                },
                "vertexgroups": {
                    "leftvertices": ["0","a"],
                    "rightvertices": ["ac", "b"],
                    "all": ["0", "a", "ac", "b", etc...]
                },
                "cellgroups": {
                    "stiffcells": ["a",etc...],
                    "squishycells": ["b", etc...],
                    "all_vertices": ["a","b", etc...]
                }
            }
        """
        
        return {
            "geometry": self._geometry_to_json(),
            "vertexgroups": self._vertex_groups_to_json(),
            "cellgroups": self._cell_groups_to_json(),
        }

class IntegratorSettings:
    @staticmethod
    def from_json(jobj):
        return IntegratorSettings(
            vertex_friction_gamma=jobj["vertex_friction_gamma"],
            step_dt=jobj["step_dt"],
        )
    
    def __init__(self, vertex_friction_gamma, step_dt):
        self._vertex_friction_gamma = vertex_friction_gamma
        self._step_dt = step_dt
    
    def to_json(self):
        return {
            "vertex_friction_gamma": self._vertex_friction_gamma,
            "step_dt": self._step_dt,
        }




class CellAreaForce:
    @staticmethod
    def from_json(jobj):
        return CellAreaForce(
            A0=jobj["A0"],
            kappa=jobj["kappa"],
        )
    
    def __init__(self, A0, kappa):
        self._A0 = A0
        self._kappa = kappa
    
    def to_json(self):
        return {
            "type": "area",
            "A0": self._A0,
            "kappa": self._kappa,
        }
    
class CellPerimeterForce:
    @staticmethod
    def from_json(jobj):
        return CellPerimeterForce(
            gamma=jobj["gamma"],
            lam=jobj["lambda"],
        )
    
    def __init__(self, gamma, lam):
        self._gamma = gamma
        self._lambda = lam
    
    def to_json(self):
        return {
            "type": "perimeter",
            "gamma": self._gamma,
            "lambda": self._lambda,
        }

class CellForce:
    @staticmethod
    def from_json(jobj):
        if jobj["type"] == "area":
            return CellAreaForce.from_json(jobj)
        elif jobj["type"] == "perimeter":
            return CellPerimeterForce.from_json(jobj)
        else:
            raise ValueError("Unknown cell force type '{}'".format(jobj["type"]))

class CellGroupForces:
    @staticmethod
    def from_json(jobj):
        c_groups_dict = jobj["cellgroups"]
        
        # cell_groups
        parsed_groups = {}
        for cell_group_name, forces_list in c_groups_dict.items():
            # forces_list = c_groups[cell_group_name]
            
            parsed_forces = []
            for force_obj in forces_list:
                cell_force_config = CellForce.from_json(force_obj)
                parsed_forces.append(cell_force_config)
            
            parsed_groups[cell_group_name] = parsed_forces
        
        return CellGroupForces(group_forces=parsed_groups)
    
    def __init__(self, group_forces):
        self._force_lists_by_groupname = group_forces
    
    def to_json(self):
        
        groups_dict = {}
        
        for grp_name, forces_for_group in self._force_lists_by_groupname.items():
            groups_dict[grp_name] = [cell_force.to_json() for cell_force in forces_for_group]
        
        return {
            'cellgroups': groups_dict,
        }
        
        
class ConstantVertexForce:
    @staticmethod
    def from_json(jobj):
        return ConstantVertexForce(
            f_x=jobj["f_x"],
            f_y=jobj["f_y"],
        )
    
    def __init__(self, f_x, f_y):
        self._f_x = f_x
        self._f_y = f_y
    
    def to_json(self):
        return {
            "type": "constant_force",
            "f_x": self._f_x,
            "f_y": self._f_y,
        }
    
class VertexForce:
    @staticmethod
    def from_json(jobj):
        if jobj["type"] == "constant_force":
            return ConstantVertexForce.from_json(jobj)
        else:
            raise ValueError("Unknown vertex force type '{}'".format(jobj["type"]))

class VertexGroupForces:
    def from_json(jobj):
        v_groups_dict = jobj["vertexgroups"]
        
        parsed_groups = {}
        for vtx_group_name, force_list in v_groups_dict.items():
            
            parsed_forces = []
            for force_obj in forces_list:
                vtx_force_config = VertexForce.from_json(force_obj)
                parsed_forces.append(vtx_force_config)
            
            parsed_groups[vtx_group_name] = parsed_forces
        
        return VertexGroupForces(vertex_forces=parsed_groups)
    
    def __init__(self, vertex_group_forces):
        self._forces_lists_by_groupname = vertex_group_forces
    
    def to_json(self):
        groups_dict = {}
        
        for grp_name, forces_for_group in self._forces_lists_by_groupname.items():
            groups_dict[grp_name] = [vtx_force.to_json() for vtx_force in forces_for_group]
        
        return {
            "vertexgroups"; groups_dict,
        }
            

class AllForces:
    @staticmethod
    def from_json(jobj):
        return AllForces(
            cell_forces=CellGroupForces.from_json(jobj["cell_forces"]),
            vertex_forces=VertexGroupForces.from_json(jobj["vertex_forces"]),
        )
    
    def __init__(self, cell_forces, vertex_forces):
        self._cell_forces = cell_forces
        self._vertex_forces = vertex_forces
    
    def to_json(self):
        return {
            "cell_forces": self._cell_forces.to_json(),
            "vertex_forces": self._vertex_forces.to_json(),
        }

class SimulationSettings:
    @staticmethod
    def from_json(jobj):
        return SimulationSettings(
            integrator_settings=IntegratorSettings.from_json(jobj["integrator"]),
            force_settings=AllForces.from_json(jobj["forces"]),
        )
    
    def __init__(self, integrator_settings, force_settings):
        self._integrator_settings = integrator_settings
        self._force_settings = force_settings
    
    def to_json(self):
        return {
            "integrator_settings": self._integrator_settings.to_json(self),
            "forces": self._force_settings.to_json(self),
        }


class VMState:
    @staticmethod
    def from_json(jobj):
        return VMState(
            tiss_topology=TissueTopology.from_json(jobj["topology"])
            current_state=TissueState.from_json(jobj["current_state"]),
            sim_settings=SimulationSettings.from_json(jobj["simsettings"])
        )
        
    def __init__(self, tiss_topology, current_state, sim_settings):
        self._tissue_topology = tiss_topology
        self._current_state = current_state
        self._simulation_settings = sim_settings
    
    def to_json(self):
        return {
            "topology": self._tissue_topology.to_json(),
            "current_state": self._current_state.to_json(),
            "simsettings": self._simulation_settings.to_json(),
        }
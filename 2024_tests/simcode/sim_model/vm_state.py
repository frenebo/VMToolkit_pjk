

# class CellTopology:
#     @staticmethod
#     def from_josn(?)

class CellTopology:
    @staticmethod
    def from_json(jobj):
        vertex_ids = jobj["vertex_ids"]
        return CellTopology(
            vertices=vertex_ids,
            is_outer=jobj["is_outer"]
        )
        
    def __init__(self, vertices, is_outer):#, is_boundary):
        self._vertex_ids=vertices
        self._is_outer=is_outer
        # self._is_boundary = is_boundary
    
    def vertex_ids(self):
        return self._vertex_ids
    
    def is_outer(self):
        return self._is_outer
    
    
    def to_json(self):
        return {
            "vertex_ids": self._vertex_ids,
            "is_outer": self._is_outer,
            # "is_boundary": self._is_boundary,
        }

class VertexTopology:
    @staticmethod
    def from_json(jobj):
        return VertexTopology(
            is_boundary=jobj["is_boundary"],
            # num_faces=jobj["num_faces"],
        )
    
    def __init__(self, is_boundary):
        self._is_boundary = is_boundary
        # pass
    
    def is_boundary(self):
        return self._is_boundary
    
    # def num_faces(self):
    #     return self._num_faces
    
    def to_json(self):
        return {
            "is_boundary": self._is_boundary,
            # "num_faces": self._num_faces
        }

class TissueTopology:
    @staticmethod
    def from_json(jobj):
        cell_topologies_parsed = {}
        
        for cell_id, cell_top_json in jobj["cells"].items():
            cell_topologies_parsed[cell_id] = CellTopology.from_json(cell_top_json)
        
        vertex_topologies_parsed = {}
        # if "vertices": not in Tissue
        for vtx_id, vtx_top_json in jobj["vertices"].items():
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
    
    def cells(self):
        return self._cell_topologies
    
    def vertices(self):
        return self._vertex_topologies
        
    def __init__(self, cell_topologies, vertex_topologies):
        self._cell_topologies = cell_topologies
        self._vertex_topologies = vertex_topologies
    
    
    def to_json(self):
        json_cell_tops = {}
        for cell_id, cell_top in self._cell_topologies.items():
            json_cell_tops[cell_id] = cell_top.to_json()
        
        json_vtx_tops = {}
        for vtx_id, vtx_top in self._vertex_topologies.items():
            json_vtx_tops[vtx_id] = vtx_top.to_json()
        return {
            "cells": json_cell_tops,
            "vertices": json_vtx_tops,
        }
        

class VertexGeometry:
    @staticmethod
    def from_json(jobj):
        return VertexGeometry(
            x=jobj["x"],
            y=jobj["y"],
        )
    
    def __init__(self, x, y):
        self._x = x
        self._y = y
    
    def x(self):
        return self._x
    
    def y(self):
        return self._y
    
    def to_json(self):
        return {
            "x": self._x,
            "y": self._y,
        }


class TissGeometry:
    @staticmethod
    def from_json(jobj):
        vtx_geoms = {}
        
        for vtx_id, vtx_geo_jobj in jobj["vertices"].items():
            vtx_geoms[vtx_id] = VertexGeometry.from_json(vtx_geo_jobj)
        
        
        return TissGeometry(
            vertex_geometries=vtx_geoms
        )
    
    def __init__(self, vertex_geometries):
        self._vertex_geometries = vertex_geometries
    
    def vertices(self):
        return self._vertex_geometries
    
    def to_json(self):
        vertices_jobjdict = {}
        for vtx_id, vtx_geom in self._vertex_geometries.items():
            # self.
            vertices_jobjdict[vtx_id] = vtx_geom.to_json()
            # self.vtx
        return {
            "vertices": vertices_jobjdict,
        }

class CellGroup:
    @staticmethod
    def from_json(jobj):
        return CellGroup(
            member_cells=jobj["cellids"]
        )
    
    def __init__(self, member_cells):
        self._cell_ids = member_cells
    
    def cell_ids(self):
        return self._cell_ids
    
    def to_json(self):
        return {
            "cellids": self._cell_ids
        }

class VertexGroup:
    @staticmethod
    def from_json(jobj):
        return VertexGroup(
            vertex_ids=jobj["vtxids"],
        )
    
    def __init__(self, vertex_ids):
        self._vertex_ids = vertex_ids
    
    def vertex_ids(self):
        return self._vertex_ids
    
    def to_json(self):
        return {
            "vtxids": self._vertex_ids
        }


class TissueState:
    @staticmethod
    def from_json(jobj):
        # vertex_geoms = {}
        
        # for vid, vertex_geometry_jobj in jobj["geometry"]["vertices"]:
        
        vertex_groups = {}
        for grpid, vtx_group_jobj in jobj["vertexgroups"].items():
            vertex_groups[grpid] = VertexGroup.from_json(vtx_group_jobj)
        
        cell_groups = {}
        for grpid, cell_group_jobj in jobj["cellgroups"].items():
            cell_groups[grpid] = CellGroup.from_json(cell_group_jobj)
        return TissueState(
            geometry=TissGeometry.from_json(jobj["geometry"]),
            vertex_groups=vertex_groups,
            cell_groups=cell_groups,
        )
        
    def __init__(self, geometry, vertex_groups, cell_groups):
        self._geometry = geometry
        self._vertex_groups = vertex_groups
        self._cell_groups = cell_groups
    
    def vertex_groups(self):
        return self._vertex_groups
        
    
    def cell_groups(self):
        return self._cell_groups
    
    def geometry(self):
        return self._geometry
    
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
        
        vtx_groups_json = {}
        for grpid, vtx_grp in self._vertex_groups.items():
            vtx_groups_json[grpid] = vtx_grp.to_json()
        
        cell_groups_json = {}
        for grpid, cell_grp in self._cell_groups.items():
            cell_groups_json[grpid] = cell_grp.to_json()
        
        return {
            "geometry": self._geometry.to_json(),
            "vertexgroups": vtx_groups_json,
            "cellgroups": cell_groups_json,
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
    
    def step_dt(self):
        return self._step_dt
    
    def vtx_friction_gamma(self):
        return self._vertex_friction_gamma
    
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
    
    def A0(self):
        return self._A0
    
    def kappa(self):
        return self._kappa
    
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
    
    def gamma(self):
        return self._gamma
    
    def lam(self):
        return self._lambda
    
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
        # c_groups_dict = jobj["cellgroups"]
        
        # cell_groups
        parsed_groups = {}
        for cell_group_name, forces_list in jobj.items():
            # forces_list = c_groups[cell_group_name]
            
            parsed_forces = []
            for force_obj in forces_list:
                cell_force_config = CellForce.from_json(force_obj)
                parsed_forces.append(cell_force_config)
            
            parsed_groups[cell_group_name] = parsed_forces
        
        return CellGroupForces(group_forces=parsed_groups)
    
    def __init__(self, group_forces):
        self._force_lists_by_groupname = group_forces

    def force_lists_by_group(self):
        return self._force_lists_by_groupname
    
    def to_json(self):
        
        groups_dict = {}
        
        for grp_name, forces_for_group in self._force_lists_by_groupname.items():
            groups_dict[grp_name] = [cell_force.to_json() for cell_force in forces_for_group]
        
        return groups_dict
        
        
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
    
    def f_x(self):
        return self._f_x
    
    def f_y(self):
        return self._f_y
        
    
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
        
        parsed_groups = {}
        for vtx_group_name, forces_list in jobj.items():
            
            parsed_forces = []
            for force_obj in forces_list:
                vtx_force_config = VertexForce.from_json(force_obj)
                parsed_forces.append(vtx_force_config)
            
            parsed_groups[vtx_group_name] = parsed_forces
        
        return VertexGroupForces(vertex_group_forces=parsed_groups)
    
    def __init__(self, vertex_group_forces):
        self._force_lists_by_groupname = vertex_group_forces
    
    def force_lists_by_group(self):
        return self._force_lists_by_groupname
    
    def to_json(self):
        groups_dict = {}
        
        for grp_name, forces_for_group in self._force_lists_by_groupname.items():
            groups_dict[grp_name] = [vtx_force.to_json() for vtx_force in forces_for_group]
        
        return groups_dict
            

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
    
    def cell_group_forces(self):
        return self._cell_forces
    
    def vertex_group_forces(self):
        return self._vertex_forces

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
    
    def integrator_settings(self):
        return self._integrator_settings
    
    def force_settings(self):
        return self._force_settings
        
        
    def to_json(self):
        return {
            "integrator": self._integrator_settings.to_json(),
            "forces": self._force_settings.to_json(),
        }


class VMState:
    @staticmethod
    def from_json(jobj):
        return VMState(
            tiss_topology=TissueTopology.from_json(jobj["topology"]),
            current_state=TissueState.from_json(jobj["current_state"]),
            sim_settings=SimulationSettings.from_json(jobj["simsettings"]),
        )
        
    def __init__(self, tiss_topology, current_state, sim_settings):
        self._tissue_topology = tiss_topology
        self._current_state = current_state
        self._simulation_settings = sim_settings
    
    def topology(self):
        return self._tissue_topology
    
    def current_state(self):
        return self._current_state
    
    def sim_settings(self):
        return self._simulation_settings
    
    def to_json(self):
        return {
            "topology": self._tissue_topology.to_json(),
            "current_state": self._current_state.to_json(),
            "simsettings": self._simulation_settings.to_json(),
        }
from .force_specs import TissueForce
import copy

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
            member_cells=jobj["cell_ids"],
            force_ids=jobj["force_ids"],
        )
    
    def __init__(self, member_cells, force_ids):
        self._cell_ids = member_cells
        self._force_ids = force_ids
    
    def cell_ids(self):
        return self._cell_ids
    
    def force_ids(self):
        return self._force_ids
    
    def set_force_ids(self, new_force_ids):
        assert isinstance(new_force_ids, list)
        for force_id in new_force_ids:
            assert type(force_id) == str
        
        self._force_ids = list(new_force_ids)
    
    def to_json(self):
        return {
            "cell_ids": list(self._cell_ids),
            "force_ids": list(self._force_ids),
        }

class VertexGroup:
    @staticmethod
    def from_json(jobj):
        return VertexGroup(
            vertex_ids=jobj["vtxids"],
            force_ids=jobj["force_ids"],
        )
    
    def __init__(self, vertex_ids, force_ids):
        self._vertex_ids = vertex_ids
        self._force_ids = force_ids
    
    def vertex_ids(self):
        return self._vertex_ids
    
    
    def force_ids(self):
        return self._force_ids
    
    def set_force_ids(self, new_force_ids):
        assert isinstance(new_force_ids, list)
        for force_id in new_force_ids:
            assert type(force_id) == str
        
        self._force_ids = list(new_force_ids)
    
    def to_json(self):
        return {
            "vtxids": list(self._vertex_ids),
            "force_ids": list(self._force_ids),
        }

class TissueState:
    @staticmethod
    def from_json(jobj):
        vertex_groups = {}
        for grpid, vtx_group_jobj in jobj["vertexgroups"].items():
            vertex_groups[grpid] = VertexGroup.from_json(vtx_group_jobj)
        
        cell_groups = {}
        for grpid, cell_group_jobj in jobj["cellgroups"].items():
            cell_groups[grpid] = CellGroup.from_json(cell_group_jobj)
        
        forces = {}
        for forceid, force_obj in jobj["forces"].items():
            forces[forceid] = TissueForce.from_json(force_obj)
        
        # if jobj["current_experienced_forces"] is None:
        #     current_experienced_forces = None
        # else:
        #     current_experienced_forces = CurrentExperiencedForces.from_json(
        #         jobj["current_experienced_forces"]
        #     )
            
        return TissueState(
            geometry=TissGeometry.from_json(jobj["geometry"]),
            vertex_groups=vertex_groups,
            cell_groups=cell_groups,
            forces=forces,
            # current_experienced_forces=current_experienced_forces,
        )
        
    def __init__(
        self,
        geometry,
        vertex_groups,
        cell_groups,
        forces,
        # current_experienced_forces,
        ):
        self._geometry = geometry
        self._vertex_groups = vertex_groups
        self._cell_groups = cell_groups
        self._forces = forces
        # self._current_experienced_forces = current_experienced_forces
    
    def vertex_groups(self):
        return self._vertex_groups
        
    
    def cell_groups(self):
        return self._cell_groups
    
    def geometry(self):
        return self._geometry

    def forces(self):
        return self._forces
    
    
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
            
        forces_json = {}
        for forceid, force in self._forces.items():
            forces_json[forceid] = force.to_json()
        
        return {
            "geometry": self._geometry.to_json(),
            "vertexgroups": vtx_groups_json,
            "cellgroups": cell_groups_json,
            "forces": forces_json,
        }
        

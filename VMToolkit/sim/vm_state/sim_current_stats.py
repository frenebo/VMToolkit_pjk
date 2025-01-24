import copy

class ForceStatsForVertices:
    @classmethod
    def from_json(cls, jobj):
        assert isinstance(jobj, dict)
        for vid, vforce in jobj["force_experienced_by_vertex_id"].items():
            assert isinstance(vforce, list)
            assert len(vforce) == 2, "should be a list of two elements, x and y force"
            for el in vforce:
                assert isinstance(el, int) or isinstance(el, float)
        return ForceStatsForVertices(
            force_experienced_by_vertex_id=jobj["force_experienced_by_vertex_id"]
        )
    
    def __init__(self, force_experienced_by_vertex_id):
        self._force_experienced_by_vertex_id = force_experienced_by_vertex_id
    
    def force_experienced_by_vertex_id(self):
        return self._force_experienced_by_vertex_id
    
    def to_json(self):
        return {
            "force_experienced_by_vertex_id": copy.deepcopy(self._force_experienced_by_vertex_id)
        }
        

class VertexAllForcesStats:
    @classmethod
    def from_json(cls, jobj):
        vtx_force_stats_by_forceid = {}
        for forceid, vertices_force_stats_jobj in jobj.items():
            vtx_force_stats_by_forceid[forceid] = ForceStatsForVertices.from_json(vertices_force_stats_jobj)
        
        return VertexAllForcesStats(
            vtx_force_stats_by_forceid,
        )
    
    def __init__(
        self,
        vertex_force_stats_by_forceid,
    ):
        assert isinstance(vertex_force_stats_by_forceid, dict)
        for forceid, vtx_forcestats in vertex_force_stats_by_forceid.items():
            assert isinstance(vtx_forcestats, ForceStatsForVertices)
        
        self._vertex_force_stats_by_forceid = vertex_force_stats_by_forceid
    
    def vertex_force_stats_by_forceid(self):
        return self._vertex_force_stats_by_forceid
    
    def to_json(self):
        vertex_force_stats_by_forceid_jobjs = {}
        
        for forceid, vertices_force_stats in self._vertex_force_stats_by_forceid.items():
            vertex_force_stats_by_forceid_jobjs[forceid] = vertices_force_stats.to_json()
            
        return vertex_force_stats_by_forceid_jobjs
        

class CurrentForcesStats:
    @classmethod
    def from_json(cls, jobj):
        return CurrentForcesStats(
            current_vertex_forces=VertexAllForcesStats.from_json(
                jobj["current_vertex_forces"]
            ),
        )
    
    def __init__(self, current_vertex_forces):
        self._current_vertex_forces = current_vertex_forces
    
    def current_vertex_forces(self):
        return self._current_vertex_forces
    
    def to_json(self):
        return {
            "current_vertex_forces": self._current_vertex_forces.to_json()
        }

class SimCurrentStats:
    @classmethod
    def from_json(cls, jobj):
        return SimCurrentStats(
            forces_stats=CurrentForcesStats.from_json(jobj["forces_stats"])
        )
    
    def __init__(self, forces_stats):
        assert isinstance(forces_stats, CurrentForcesStats)
        
        self._forces_stats = forces_stats
    
    def forces_stats(self):
        return self._forces_stats
    
    def to_json(self):
        return {
            "forces_stats": self._forces_stats.to_json()
        }
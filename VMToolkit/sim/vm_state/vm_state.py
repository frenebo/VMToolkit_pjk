
from .simulation_settings import SimulationSettings
from .tissue_topology import TissueTopology
from .tissue_state import TissueState
from .sim_current_stats import SimCurrentStats
import copy

class VMState:
    @staticmethod
    def from_json(jobj):
        if jobj["sim_current_stats"] is None:
            sim_current_stats = None
        else:
            sim_current_stats = SimCurrentStats.from_json(jobj["sim_current_stats"])
        return VMState(
            tiss_topology=TissueTopology.from_json(jobj["topology"]),
            current_tissue_state=TissueState.from_json(jobj["current_tissue_state"]),
            sim_settings=SimulationSettings.from_json(jobj["simsettings"]),
            sim_current_stats=sim_current_stats
        )
        
    def __init__(self, tiss_topology, current_tissue_state, sim_settings, sim_current_stats):
        self._tissue_topology = tiss_topology
        self._current_tissue_state = current_tissue_state
        self._simulation_settings = sim_settings
        self._sim_current_stats = sim_current_stats
    
    def topology(self):
        return self._tissue_topology
    
    def tissue_state(self):
        return self._current_tissue_state
    
    def sim_settings(self):
        return self._simulation_settings
    
    def sim_current_stats(self):
        return self._sim_current_stats
    
    def set_sim_current_stats(self, sim_current_stats):
        assert isinstance(sim_current_stats, SimCurrentStats) or sim_current_stats is None
        
        self._sim_current_stats = sim_current_stats
    
    def to_json(self):
        if self._sim_current_stats is None:
            sim_current_stats_jobj = None
        else:
            sim_current_stats_jobj = self._sim_current_stats.to_json()
        return {
            "topology": self._tissue_topology.to_json(),
            "current_tissue_state": self._current_tissue_state.to_json(),
            "simsettings": self._simulation_settings.to_json(),
            "sim_current_stats": sim_current_stats_jobj,
        }


from .simulation_settings import SimulationSettings
from .tissue_topology import TissueTopology
from .tissue_state import TissueState
import copy

class VMState:
    @staticmethod
    def from_json(jobj):
        return VMState(
            tiss_topology=TissueTopology.from_json(jobj["topology"]),
            current_tissue_state=TissueState.from_json(jobj["current_tissue_state"]),
            sim_settings=SimulationSettings.from_json(jobj["simsettings"]),
        )
        
    def __init__(self, tiss_topology, current_tissue_state, sim_settings):
        self._tissue_topology = tiss_topology
        self._current_tissue_state = current_tissue_state
        self._simulation_settings = sim_settings
    
    def topology(self):
        return self._tissue_topology
    
    def tissue_state(self):
        return self._current_tissue_state
    
    def sim_settings(self):
        return self._simulation_settings
    
    def to_json(self):
        return {
            "topology": self._tissue_topology.to_json(),
            "current_tissue_state": self._current_tissue_state.to_json(),
            "simsettings": self._simulation_settings.to_json(),
        }

import json
# from .vm_state import 
from .vm_tk_wrapper import VMToolkitWrapper
from .vm_state import VMState




# @TODO the categories and forces should be totally abstracted out of the base vertex model and wrapper
# the vertex model and the cass that controls it should only care about the forces for each individual vertex.
class SimModel:
    def __init__(self, verbose):
        self._vm_wrapper = VMToolkitWrapper(verbose=verbose)
        
        self._last_vm_state = None
    
    def load_from_json_state(self, vm_state_json)
        vm_state = VMState.from_json(vm_state_json)
        self._vm_wrapper.initialize_from_vm_state(vm_state)
    
    def configure_integrators(self,verbose=False):
        self._vm_wrapper.configure_integrators(verbose=verbose)
    
    def dump_cpp_json(self):
        return self._vm_wrapper.dump_cpp_json()
        
    def get_json_state(self):
        vertices = {}
        
        cpp_json = json.loads(self._vm_wrapper.dump_cpp_json())
    
        for vtx_index, vtx in enumerate(cpp_json["mesh"]["vertices"]):
            print(vtx)
            vtx_id = str(vtx_index)
            vertices[vtx_id] = {
                "x": vtx['r'][0],
                "y": vtx['r'][1],
            }
        
        
        # edges = {}
        cells = {}
        for cell_index, cpp_cell in enumerate(cpp_json["mesh"]["faces"]):
            cell_id = str(cell_index)
            cells[cell_id] = {
                "vertices": [str(vidx) for vidx in cpp_cell["vertices"]]
            }
        
        return {
            "vertices": vertices,
            "cells": cells,
        }

    def configure_forces(self, P0_model, gamma, kappa, verbose=False):
        self._vm_wrapper.configure_forces(P0_model, gamma, kappa, verbose=verbose)
    
    def load_json_obj(self, json_obj, verbose=False):
        self._vm_wrapper.load_json_obj(json_obj, verbose=verbose)

    
        
    def run_steps(self, n_steps, verbose=False):
        self._vm_wrapper.run_steps(n_steps, verbose=verbose)
        
        

default_config = {
    "integrator": {
        "type": "brownian",
        "temp": 0,
        "friction_gam": 1.0,
        "dt": 0.3,
    },
    
    "cell_types": {
        "default": {
            "forces": {
                "area": {
                    'kappa': 1.0,
                },
                "perimeter": {
                    'gamma': 0.15,
                    'lambda': 0.15 * 0.20 ,
                },
            },
        },
    }
}
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
    
    def load_from_json_state(self, vm_state_json, verbose=False):
        vm_state = VMState.from_json(vm_state_json)
        self._vm_wrapper.initialize_from_vm_state(vm_state, verbose=verbose)
    
    
    def dump_cpp_json(self):
        return self._vm_wrapper.dump_cpp_json()
     
        
    def run_steps(self, n_steps, verbose=False):
        self._vm_wrapper.run_steps(n_steps, verbose=verbose)
    
    # def select_g
        
        

# default_config = {
#     "integrator": {
#         "type": "brownian",
#         "temp": 0,
#         "friction_gam": 1.0,
#         "dt": 0.3,
#     },
    
#     "cell_types": {
#         "default": {
#             "forces": {
#                 "area": {
#                     'kappa': 1.0,
#                 },
#                 "perimeter": {
#                     'gamma': 0.15,
#                     'lambda': 0.15 * 0.20 ,
#                 },
#             },
#         },
#     }
# }
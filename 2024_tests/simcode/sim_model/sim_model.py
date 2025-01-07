import json
# from .vm_state import 
from .vm_wrapper.vm_tk_wrapper import VMToolkitWrapper
from .vm_state import VMState
from .model_change_requests import ModelChangeRequest
from .model_change_req_applicator import ModelChangeReqApplicator




# @TODO the categories and forces should be totally abstracted out of the base vertex model and wrapper
# the vertex model and the cass that controls it should only care about the forces for each individual vertex.
class SimModel:
    def __init__(self, verbose):
        self._vm_wrapper = VMToolkitWrapper(verbose=verbose)
    
    def load_from_json_state(self, vm_state_json, verbose=False):
        vm_state = VMState.from_json(vm_state_json)
        self._vm_wrapper.initialize_from_vm_state(vm_state, verbose=verbose)
    
    def vm_state_json(self):
        return self._vm_wrapper.vm_state_json()
        
    def run_steps(self, n_steps, verbose=False, do_time_force_computation=False):
        if do_time_force_computation:
            self._vm_wrapper.start_force_compute_timers(verbose=verbose)

        self._vm_wrapper.run_steps(n_steps, verbose=verbose)
        if do_time_force_computation:
            print("Force computation times (milliseconds):")
            print(self._vm_wrapper.get_force_compute_timers_millis(verbose=verbose))
        
    def _check_topology_changed(self):
        return self._vm_wrapper._check_topology_changed()
        
    
    def process_model_change_req(self, model_change_req):
        assert isinstance(model_change_req, ModelChangeRequest)
        
        raise NotImplementedError()
        # try:
        #     ModelChangeReqApplicator.apply_request_to_vm_state(self._last_vm_state, model_change_req)
        # except:
        #     print("Failed to apply model change rquest to current vm state.")
        #     print("Model change req: {}".format(model_change_req))
        #     print("Model change req to json: {}".format(model_change_req.to_json()))
        #     raise

import json

from .vm_wrapper import VMCppWrapper
from .vm_state import VMState




# @TODO the categories and forces should be totally abstracted out of the base vertex model and wrapper
# the vertex model and the cass that controls it should only care about the forces for each individual vertex.
class SimModel:
    def __init__(self, verbose):
        self._vm_wrapper = VMCppWrapper(verbose=verbose)
    
    def load_from_json_state(self, vm_state_json, verbose=False):
        vm_state = VMState.from_json(vm_state_json)
        self._vm_wrapper.initialize_from_vm_state(vm_state, verbose=verbose)
    
    def vm_state_json(self):
        return self._vm_wrapper.vm_state_json()
        
    def run_steps_manual_tstep(self, n_steps, verbose=False, do_time_force_computation=False):
        if do_time_force_computation:
            self._vm_wrapper.start_force_compute_timers(verbose=verbose)

        self._vm_wrapper.run_steps_manual_tstep(n_steps, verbose=verbose)
        if do_time_force_computation:
            print("Force computation times (milliseconds):")
            print(self._vm_wrapper.get_force_compute_timers_millis(verbose=verbose))
            
    def run_with_adaptive_tstep(self, time_run, verbose=False, do_time_force_computation=False):
        if do_time_force_computation:
            self._vm_wrapper.start_force_compute_timers(verbose=verbose)

        self._vm_wrapper.run_with_adaptive_tstep(time_run, verbose=verbose)
        if do_time_force_computation:
            print("Force computation times (milliseconds):")
            print(self._vm_wrapper.get_force_compute_timers_millis(verbose=verbose))
        
    def _check_topology_changed(self):
        return self._vm_wrapper._check_topology_changed()
        
    
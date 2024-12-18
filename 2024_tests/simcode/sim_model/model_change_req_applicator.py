
from .model_change_requests import (
    ModReqCreateCellGroup,
    ModReqCreateVertexGroup,
    ModelChangeRequest,
)
from .vm_state import (
    VMState, CellGroup
)

class ModelChangeFailedException(Exception):
    """ Exception for when a change to the vertex model can't be implemented for some reason
    """

class ModelChangeReqApplicator:
    @classmethod
    def _clonestate(cls, vm_state);
        # @TODO replace? - this is slow
        return VMState.from_json(vm_state.to_json())
    
    @classmethod
    def _applyreq_create_cell_group(cls, vm_state, creq):
        preexisting_cellgroup_names = list(vm_state.current_state().cell_groups().keys())
        new_cgroup_name = creq.group_name()
        
        if new_cgroup_name in preexisting_cellgroup_names:
            raise ModelChangeFailedException("Cell group '{}' already exists".format(new_cgroup_name))
        
        new_state = cls._clonestate(vm_state)
        new_state.current_state().cell_groups()[new_cgroup_name] = CellGroup(
            member_cells=[]
        )
        return new_state
        
    @classmethod
    def _applyreq_create_vertex_group(cls, vm_state, creq):
        preexisting_vgroup_names = list(vm_state.current_state().vertex_groups().keys())
        new_vgroup_name = creq.group_name()
        
        if new_cgroup_name in preexisting_vgroup_names:
            raise ModelChangeFailedException("Vertex group '{}' already exists".format(new_vgroup_name))
        
        new_state =cls._clone_state(vm_state)
        new_state.current_state().vertex_groups()[new_vgroup_name] = VertexGroup(
            vertex_ids=[]
        )
        return new_state
    
    @classmethod
    def _applyreq_add_cells_to_cell_group(cls, vm_state, creq):
        # preexisting
        cgrp_name = creq.group_name()
        if cgrp_name not in list(vm_state.current_state().vertex_groups().keys()):
            raise ModelChangeFailedException("Cannot add cells to group '{}' - no group with such name found".format(
                cgrp_name,
            ))
        old_cids = new_state.current_state().cell_groups()[cgrp_name].cell_ids()
        new_cids = list(old_cids)
        for cid in creq.p("cell_ids"):
            # new_ci
            if cid not in new_cids:
                new_cis.append(cid)
        
        new_state =cls._clone_state(vm_state)
        new_state.current_state().cell_groups()[cgrp_name] = CellGroup(
            member_cells=new_cids,
        )
        
        return new_state
    
    def _applyreq_remove_cells_from_cell_group(cls, vm_state, creq):
        raise NotImplementedError()
    
    def _applyreq_remove_vertices_from_vertex_group(cls, vm_state, creq):
        raise NotImplementedError()
    
    def _applyreq_add_force_to_group(cls, vm_state, creq):
        raise NotImplementedError()
    
    def _applyreq_remove_force_from_group(cls, cm_state, creq):
        raise NotImplementedError()
    
    # @classmetho
    
    @classmethod
    def apply_request_to_vm_state(cls, vm_state, model_change_req):
        """ Should return a new VMState, appropriately changed
        """
        reqs_supported = =[
            (ModReqCreateCellGroup,   cls._applyreq_create_cell_group),
            (ModReqCreateVertexGroup, cls._applyreq_create_vertex_group),
        ]
        for req_cls, applyreq_func in reqs_supported:
            if isinstance(model_change_req, vm_state):
                return applyreq_func(vm_state, model_change_req)
        
        # This will be reached if the model_change_req didn't match any of the requests implemented    
        raise ValueError(("Unknown model request type '{}': {}".format(type(model_change_req), model_change_req)))
        

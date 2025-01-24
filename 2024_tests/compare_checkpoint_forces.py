
import json
import argparse
from VMToolkit.sim.vm_state import (
    VMState,
    # SimulationSettings,
    # IntegratorSettings, 
    # CellAreaForce,
    # CellPerimeterForce,
    # ConstantVertexForce,
    # CellGroup,
    # VertexGroup,
    # UniformElectricForceOnCellBoundary,
    # EFieldSpecConstantPolygonRegion,
    # PolygonSpec,
    # PixelatedElectricForceOnCellBoundary,
    # PixElecForceCellParam,
    # PixelatedFieldSpec,
    # TopologySettings,
    # T1TransitionSettings,
    # IntegratorSettingsDormandPrinceRungeKutta,
)

def get_cells_affected_by_force(ckpt, forceid):
    cids_affected = []
    for cellgroup_name, cgroup in ckpt.tissue_state().cell_groups().items():
        if forceid not in cgroup.force_ids():
            continue
        # print(cgroup.cell_ids())
        # print(cgroup.force_ids())
        cids_affected += cgroup.cell_ids()
    cids_affected = list(set(cids_affected))
    
    return cids_affected
    # ckpt2_state.tissue_state().force
    
def get_vtx_forces_from_forceid(state, forceid):
    return state.sim_current_stats().forces_stats().current_vertex_forces().vertex_force_stats_by_forceid()[forceid].force_experienced_by_vertex_id() 
def get_forcetype(cstate, forceid):
    ftype = cstate.tissue_state().forces()[forceid].force_type
    print(ftype)
    return ftype
def compare_ckpt_forces(
    ckpt1_fpath,
    ckpt2_fpath,
    forceid,
):
    print("Comparing forces from '{}'".format(forceid))
    with open(ckpt1_fpath, "r") as f:
        ckpt1_state = VMState.from_json(json.load(f))
    
    with open(ckpt2_fpath, "r") as f:
        ckpt2_state = VMState.from_json(json.load(f))
    print("{ftype1} , {ftype2}".format(
        ftype1=get_forcetype(ckpt1_state, forceid),
        ftype2=get_forcetype(ckpt2_state, forceid),
    ))
    
    # print(ckpt1_state)s()[forceid]
    cellids_affected = sorted(get_cells_affected_by_force(ckpt1_state, forceid))
    if cellids_affected != sorted(get_cells_affected_by_force(ckpt2_state, forceid)):
        print(cellids_affected)
        print(get_cells_affected_by_force(ckpt2_state, forceid).sorted())
        raise ValueError("Cell ids don't match for cells affected by {forceid} in ckpt1 and ckpt2...")
    # print(cellids_affected)
    
    state1_vtx_forces = get_vtx_forces_from_forceid(ckpt1_state, forceid)
    state2_vtx_forces = get_vtx_forces_from_forceid(ckpt2_state, forceid)
    # print()
    # print(get_vtx_forces_from_forceid(ckpt2_state, forceid))
    
    for vtx_id in state1_vtx_forces:
        force1 = state1_vtx_forces[vtx_id]
        force2 = state2_vtx_forces[vtx_id]
        if force1 == [0.0,0.0] and force2 == [0.0,0.0]:
            continue
        print("vtx {vid} - ( {f1x: .8f},{f1y: .8f} ) :: ( {f2x: .8f},{f2y: .8f} )".format(
            vid=("'"+vtx_id+"'").ljust(6),
            f1x=force1[0],
            f1y=force1[1],
            f2x=force2[0],
            f2y=force2[1],
        ))
        # print()
        # print()
    # print(ckpt1_vtx_forces)
    # print(ckpt2_vtx_forces)
    


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    
    
    parser.add_argument("ckpt1")
    parser.add_argument("ckpt2")
    
    parser.add_argument("--forceid", type=str, default='left_forcing_field')
    
    args = parser.parse_args()
    
    compare_ckpt_forces(
        ckpt1_fpath=args.ckpt1,
        ckpt2_fpath=args.ckpt2,
        forceid=args.forceid,
    )
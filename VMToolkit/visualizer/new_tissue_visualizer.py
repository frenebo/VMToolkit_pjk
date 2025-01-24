from ..sim.vm_state import VMState


def visualize_simulation(tissue_checkpoint_fps):
    tissue_checkpoint_fps = sorted(tissue_checkpoint_fps)
    ckpt_states = []
    for ckpt_fp in tissue_checkpoint_fps:
        with open(ckpt_fp, r) as f:
            ck_state = VMState.from_json(json.load(f))
            ckpt_states.append(ck_state)
    
    print(ckpt_states)
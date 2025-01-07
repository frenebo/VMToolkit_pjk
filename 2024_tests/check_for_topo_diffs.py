import argparse
import json
import os

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    
    parser.add_argument("input_dir")
    
    args = parser.parse_args()
    
    
    if not os.path.exists(args.input_dir):
        raise ValueError("No such dir {}".format(args.input_dir))
    
    
    fns = sorted(os.listdir(args.input_dir))
    fns = [f for f in fns if f.endswith(".json")]
    
    all_topos = []
    
    for fn in fns:
        with open(os.path.join(args.input_dir, fn)) as f:
            vmstate = json.load(f)
        
        ckpt_topo = vmstate["topology"]
        all_topos.append(ckpt_topo)
        
    ### Check if topos change
    first_topo = all_topos[0]
    for ckpt_idx, top in enumerate(all_topos):
        for cid, ctopo in top["cells"].items():
            member_vids = ctopo['vertex_ids']
            orig_vids = first_topo['cells'][cid]["vertex_ids"]
            
            if not member_vids == orig_vids:
                print("Found difference at ckpt {}".format(ckpt_idx))
                raise Exception()
            # break
            print(".",end="")
        # break
        print("")
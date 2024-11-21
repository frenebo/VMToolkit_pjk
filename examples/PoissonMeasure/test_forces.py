import argparse

from VMToolkit.config_builder.open.make_honeycomb import create_honeycomb_json

def test_forces(scratch_dirpath):
    os.makedirs(scratch_dirpath)
    
    honeycomb_pth = os.path.join(scratch_dirpath, "start_honeycomb.json")
    
    create_honeycomb_json(honeycomb_pth)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    
    parser.add_argument("--scratch_dir", default="scratch", help="Name of scratch directory to make and create temporary files")
    args = parser.parse_args()
    
    if os.path.exists(args.scratch_dir) and not os.path.isdir(args.scratch_dir):
        raise ValueError("Cannot use {} as scratch dir, it exists but is not a directory".format(args.scratch_dir))
    
    test_forces(scratch_dirpath=args.scratch_dir)
    
    
    
    
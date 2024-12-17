import argparse
from tissue_builder.hexagonal import HexagonalCellMesh


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    
    parser.add_argument("--box_lx", type=float, default=10.0)
    parser.add_argument("--box_ly", type=float, default=10.0)
    
    args = parser.parse_args()
    cm = HexagonalCellMesh(
        1.0,
        args.box_lx,
        args.box_ly,
    )
    cm.set_all_A0(4)
    cm.set_all_P0(7.2)
    print(cm.build_vm_state())

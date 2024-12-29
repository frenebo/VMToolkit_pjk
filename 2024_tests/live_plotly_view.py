import sys
# sys.path.append("./") # go to parent dir
# sys.path.append("/Users/paulkreymborg/Documents/code/VMTutorial_pjk/2024_tests")
# print(sys.path)
import os
import pandas as pd
import plotly.express as px
import json
# import code
from simcode.plot_cells import visualize
# import code.plot_cells.visualize
import numpy as np

# code.plot_cells
# print(dir(code))


if __name__ == "__main__":

    ckpts_dir = "./scratch/"
    tiss_ckpt_fps = []
    for fn in sorted(os.listdir(ckpts_dir)):
        # print(fn)
        if fn.endswith(".json") and fn.startswith("res"):
            tiss_ckpt_fps.append(os.path.join(ckpts_dir, fn))

    fig = visualize.make_plotly_visualizer(
        tiss_ckpt_fps,
        "./scratch/initial_vm_state.json"
        # vertices_to_highlight=[0],
    )
    
    fig.show(renderer='browser')

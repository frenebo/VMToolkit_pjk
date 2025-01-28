from skbuild import setup

setup(
    name="VMToolkit",
    version="0.1",
    description="Vertex Model Toolkit",
    author="Rastko Sknepnek",
    license="MIT",
    packages=[
        'VMToolkit',
        'VMToolkit.VM',
        'VMToolkit.visualizer',
        'VMToolkit.theoretical_models',
        'VMToolkit.sim',
        'VMToolkit.sim.utils',
        'VMToolkit.sim.tissue_builder',
        'VMToolkit.sim.vm_state',
        'VMToolkit.sim.vm_wrapper',
        ],
    cmake_install_dir='VMToolkit/VM',
    python_requires=">=3.8",
)
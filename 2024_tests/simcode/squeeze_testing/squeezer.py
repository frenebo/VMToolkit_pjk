import numpy as np
# import json
import os


from ..theoretical_model.find_regular_hexagonal_rests_area import RegHexagonalModel, find_regular_hexagonal_rests_area

from ..tissue_builder.hexagonal import HexagonalCellMesh
from ..sim_model.sim_model import SimModel
from ..sim_model.vm_state import (
    VMState, SimulationSettings, IntegratorSettings, 
    CellAreaForce, CellPerimeterForce, ConstantVertexForce,  CellGroup, VertexGroup,
    ElectricForceOnCellBoundary, EFieldSpecConstantPolygonRegion, PolygonSpec,
    
)
from ..sim_model.box_selector import BoxSelector

class SqueezeSimAnalyzer:
    
    def __init__(
        self,
        tissue_cell_properties,
        env_setup_properties,
        simulation_settings,
        
        
    ):
        pass

class SqueezerSim:
    
    @staticmethod
    def _build_squeezer_model(
        tissue_cell_properties,
        env_setup_properties,
        simulation_settings
        ):
        A0_model = tissue_cell_properties["A0"]
        P0_model = tissue_cell_properties["P0"]
        
        gamma_force_coeff = tissue_cell_properties["gamma"]
        kappa_force_coeff = tissue_cell_properties["kappa"]
        
        box_lx = env_setup_properties["box_lx"]
        box_ly = env_setup_properties["box_ly"]
        
        sim_friction_gamma = simulation_settings["friction_gamma"]
        sim_step_dt = simulation_settings["step_dt"]
        
        forcing_field_strength = env_setup_properties["forcing_field_strength"]
        field_zone_rel_width = env_setup_properties["forcing_field_relative_width"]
        
        
        analytical_predictions = RegHexagonalModel().find_elastic_props_of_hexagon(
            A_0_num=A0_model,
            P_0_num=P0_model,
            K_num=kappa_force_coeff,
            gamma_num=gamma_force_coeff,
        )
        
        rest_side_length = analytical_predictions["rest_side_length"]
        theoretical_rest_width = rest_side_length*2
        theoretical_rest_height = rest_side_length*2*np.sqrt(3)/2
        
        print("Theoretical rest width, height: {}, {}".format(theoretical_rest_width, theoretical_rest_height))
        
        
        cm = HexagonalCellMesh(
            side_length=rest_side_length,
            box_lx=box_lx,
            box_ly=box_ly,
        )
        
        tiss_topology, tiss_init_state = cm.build_vm_state()
        
        tiss_init_state.forces()["perim_f"] = CellPerimeterForce(
            gamma=gamma_force_coeff,
            lam=P0_model*gamma_force_coeff,
        )
        tiss_init_state.forces()["area_f"] = CellAreaForce(
            A0=A0_model,
            kappa=kappa_force_coeff,
        )
        
        # field_size_multiplier = 1.0
        
        tiss_init_state.forces()["left_forcing_field"] = make_forcing_field_rectangular(
            xmin= -box_lx,
            xmax= - box_lx * (0.5 - field_zone_rel_width),
            ymin=-box_ly,
            ymax=box_ly,
            field_x=forcing_field_strength,
            field_y=0.0,
        )
        tiss_init_state.forces()["right_forcing_field"] = make_forcing_field_rectangular(
            xmin=box_lx * (0.5 - field_zone_rel_width),
            xmax=box_lx,
            ymin=-box_ly,
            ymax=box_ly,
            field_x=-forcing_field_strength,
            field_y=0.0,
        )
        
        tiss_init_state.cell_groups()["normal"].force_ids().extend([
            "perim_f",
            "area_f",
            
            "left_forcing_field",
            "right_forcing_field",
        ])
        
        vm_initial_state = VMState(
            tiss_topology=tiss_topology,
            current_tissue_state=tiss_init_state,
            sim_settings=SimulationSettings(
                integrator_settings=IntegratorSettings(
                    vertex_friction_gamma=sim_friction_gamma,
                    step_dt=sim_step_dt,
                ),
            )
        )
        
        print("MODEL PARAMS")
        print("A0={A0}  P0={P0}  kappa={kappa} gamma={gamma}".format(A0=A0_model,P0=P0_model,kappa=kappa_force_coeff,gamma=gamma_force_coeff))
        print(analytical_predictions)
        print("theoretical rest width={}, height={}".format(theoretical_rest_width, theoretical_rest_height))
            
        return vm_initial_state


    def __init__(
        self,
        tissue_cell_properties,
        env_setup_properties,
        simulation_settings,
    ):
        self._tissue_cell_properties = tissue_cell_properties
        self._env_setup_properties = env_setup_properties
        self._simulation_settings = simulation_settings
        
        # self._vm_state = 
        init_vm_state = SqueezerSim._build_squeezer_model(
            tissue_cell_properties,
            env_setup_properties,
            simulation_settings,
        )
        
        self.sim_model = SimModel(verbose=False)
        self.sim_model.load_from_json_state(init_vm_state.to_json())
        
        self._finished_sim = False
    
    def run_sim(self):
        step_size = self._simulation_settings["checkpoint_step_size"]
        n_checkpoints = self._simulation_settings["n_checkpoints"]
        
        vm_states_over_time = []
        
        vm_states_over_time.append(self.sim_model.vm_state_json())
    
        for i in range(N_checkpoints):
            sim_model.run_steps(step_size)
            state_json_obj = sim_model.vm_state_json()
        
        self._finished_sim = True
    
    def analyze_sim(self):
        if not self._finished_sim:
            raise Exception("Need to run sim before analyzing results!")
        
        
        

def make_forcing_field_rectangular(
    xmin,
    xmax,
    ymin,
    ymax,
    field_x,
    field_y,
    ):
    return ElectricForceOnCellBoundary(electric_field_spec=(
        EFieldSpecConstantPolygonRegion(
            E_x=field_x,
            E_y=field_y,
            zone_bounds=PolygonSpec(
                polygon_vertices=[
                    (xmin,ymin),
                    (xmax,ymin),
                    (xmax,ymax),
                    (xmin,ymax),
                ]
            )
        )
    ))



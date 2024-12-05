
from hexagonal_test.analytic_tools.find_hexagon_rest_area import HexagonalModel, find_hexagon_rest_area

from VMToolkit.config_builder.open.make_honeycomb import create_honeycomb_json
from VMToolkit.config_builder.open.honeycomb_lattice import HoneycombLattice


if __name__ == "__main__":
    hex_model = HexagonalModel()
    
    A0_model = 20.0
    P0_model = 1.0
    res = hex_model.find_rest_size_of_hexagon(
        A_0=A0_model,
        P_0=P0_model,
        K=1,
        Gamma=1,
    )
    rest_side_length = res["rest_side_length"]
    print(res)
    
    h = HoneycombLattice(
        lx=25.0,
        ly=25.0,
        a=rest_side_length,
    )
    h.build()
    h.minmax()
    h.set_energy_params(A0=A0_model,P0=P0_model)    
    h.json_out("scratch/example.json")
    
    
    # create_honeycomb_json(
    #     honeycomb_pth,
    #     honeycomb_a_len=side_len_rest,
    # )

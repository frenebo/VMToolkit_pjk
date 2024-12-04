
from hexagonal_test.analytic_tools.find_hexagon_rest_area import HexagonalModel, find_hexagon_rest_area

from VMToolkit.config_builder.open.make_honeycomb import create_honeycomb_json
from VMToolkit.config_builder.open.honeycomb_lattice import HoneycombLattice


if __name__ == "__main__":
    hex_model = HexagonalModel()
    res = hex_model.find_rest_size_of_hexagon(
        20,
        1,
        1,
        1,
    )
    rest_side_length = res["rest_side_length"]
    print(res)
    
    honeycomb = HoneycombLattice(
        lx=25.0,
        ly=25.0,
        a=rest_side_length,
    )
    
    honeycomb.json_out("scratch/example.json")
    
    
    # create_honeycomb_json(
    #     honeycomb_pth,
    #     honeycomb_a_len=side_len_rest,
    # )
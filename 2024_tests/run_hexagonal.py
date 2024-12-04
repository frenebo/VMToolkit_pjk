
from hexagonal_test.analytic_tools.find_hexagon_rest_area import HexagonalModel, find_hexagon_rest_area

from VMToolkit.config_builder.open.make_honeycomb import create_honeycomb_json


if __name__ == "__main__":
    hex_model = HexagonalModel()
    res = hex_model.find_rest_size_of_hexagon(
        20,
        1,
        1,
        1,
    )
    
    print(res)
    
    # create_honeycomb_json(
    #     honeycomb_pth,
    #     honeycomb_a_len=side_len_rest,
    # )
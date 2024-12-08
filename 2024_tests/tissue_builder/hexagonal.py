import math
import numpy as np


class HexagonalCellMesh:
    
    # @classmethod
    # def get_cell_row_xpositions(side_length, box_lx)
    
    # def generate_vertex_row(n_cells_horizontal)
    @classmethod
    def generate_vertex_row_untrimmed(cls, side_length, row_is_odd, num_cells_wide):
        # side
        if num_cells_wide % 2 != 1:
            raise ValueError("The tissue should be an odd number of cells wide.")
        
        vertex_x_positions = []
        
        # Two possibilities
        #########
        # Odd rows look like this - using parentheses to show middle position.:
        #    X------X--X--()--X--X------X
        #       0         1         2
        # Even rows look like this:
        #      X--X------X()X------X--X
        #       0         1         2
        
        #########
        # Odd rows look like this - using parentheses to show middle position.:
        #      X--X------X--X--()--X--X------X--X
        #       0         1         2         3
        # Even rows look like this:
        #    X------X--X------X()X------X--X------X
        #       0         1         2         3
        
        # Look at every other 'cell position', make two vertices for it.
        # range(0, num_cells_wide, 2)
        
        cell_indices_to_build_vtxs_from = list(range(0, num_cells_wide, 2))
        first_build_cell_is_leftmost = len(cell_indices_to_build_vtxs_from) % 2 == 1
        
        build_cell_vertices_wide = not (row_is_odd ^ first_build_cell_is_leftmost)
        
        vertex_x_positions = []
        
        # vertices_close_to_cells = not ( (len(cell_indices_to_build_vtxs_from) % 2 == 1) ^ row_is_odd)
        for cell_idx in cell_indices_to_build_vtxs_from:
            cell_xpos = (cell_idx - (num_cells_wide // 2 + 1)) * side_length
            # print()
            if build_cell_vertices_wide:
                left_vtx_x = cell_xpos - side_length
                right_vtx_x = cell_xpos + side_length
            else:
                left_vtx_x = cell_xpos - side_length/2
                right_vtx_x = cell_xpos + side_length/2
            
            vertex_x_positions.append(left_vtx_x)
            vertex_x_positions.append(right_vtx_x)
        
        return vertex_x_positions
        
    
    @classmethod
    def build_cell_centers(cls, side_length, box_lx, box_ly):
        hex_h = side_length*math.sqrt(3)
        hex_w = side_length*2
        print("hex_w={} hex_h={}".format(hex_w,hex_h))
        
        # The 'first' row has the height of one hexagon - each additional one adds half a hexagon height.
        # so...
        # tissue_height = (num_rows/2 + 1/2)*hexagon_height
        num_cell_rows = math.floor(2*((box_ly / hex_h) - 1/2))
        # Round down to an odd number - things will be easier if we have one cell at the origin, and build
        # outwards. otherwise, we might end up with asymettrical tissue and stuff
        num_cell_rows = ((num_cell_rows - 1)//2)*2 + 1
        
        
        # num_cell_columns = 
        #### Finding number of cell columns
        # The first cell has width hex_w. Each additional hexagon adds width 1.5*hex_w.
        # So for n_cells, total width is = hex_w + (n_cells - 1)*1.5*hex_w
        #                      width_tot = hex_w*((3/2)*n_cells - 1/2)
        #                      width_tot / (hex_w*3/2) +1/2 = n_cells
        # So finding n_cells:
        # print("Box lx={}   box_lx/(hex_w*(3/2))={}".format(box_lx, box_lx/(hex_w*3/2)))
        num_cell_columns = math.floor((box_lx / (hex_w*(3/2))) + (1/2))
        
        num_cell_columns = ( (num_cell_columns - 1)//2 )*2 + 1
        
        print("Number of cell columns: {}".format(num_cell_columns))
        print("Number of cell rows: {}".format(num_cell_rows))
        # print("Number of cell columns: {}".format(num_cell_columns))
        
        # The middle row is at origin, and the others are spaced apart from it
        cell_row_ypositions = (hex_h/2)*(np.arange(num_cell_rows) - num_cell_rows//2)
        
        # print("Cell row y positions: {}".format(cell_row_ypositions))
        
        # vtx_row_
        vtx_rows = []
        for vtx_row_idx in range(num_cell_rows + 2):
            # The middle cell row is "odd" - has cell in it's center.
            # So, the row of vertices that are along the middle of that row are 'odd'.
            #  - this middle row of vertices has index  n_cell_rows//2 + 1
            # Since vertex rows alternate being 'even' or 'odd',
            # a row with index i is odd, if (i - (n_cell_rows//2 + 1)) % 2 == 0
            # print(vtx_row_idx)
            vtx_row_is_odd = (vtx_row_idx - (num_cell_rows//2 + 1)) % 2 == 0
            
            
            # print(vtx_row_is_odd)
            vtx_rows.append(cls.generate_vertex_row_untrimmed(
                side_length=side_length,
                row_is_odd=vtx_row_is_odd,
                num_cells_wide=num_cell_columns,
            ))
        print(vtx_rows)
        
        # for row_idx in range(num_rows):
        #     y = row_ypositions
            
            # row_cent
        
        
        
        
    def __init__(
        self,
        side_length,
        box_lx,
        box_ly,
    ):
        # raise NotImplementedError()
        
        # center_cell = 
        # cent
        # centself.
        self.side_length = side_length
        self.box_lx = box_lx
        self.box_ly = box_ly
        self.build_cells()
    
    def build_cells(self):
        centers = self.build_cell_centers(self.side_length, self.box_lx, self.box_ly)    
        # raise Exception("Have not excuded the corners for even rows that are first or lat row")
        
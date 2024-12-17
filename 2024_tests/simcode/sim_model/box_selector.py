
class BoxSelector:
    @staticmethod
    def get_cells_in_box(vm_state, x1,x2,y1,y2):
        assert x2 >= x1
        assert y2 >= y2
        
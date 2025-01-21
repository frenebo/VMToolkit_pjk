import copy

class PolygonSpec:
    @staticmethod
    def from_json(jobj):
        polygon_vertices = jobj["polygon_vertices"]
        assert isinstance(polygon_vertices, list)
        for v in polygon_vertices:
            assert isinstance(v, (list,tuple))
            assert len(v) == 2
            for coord in v:
                assert isinstance(coord, (int, float, complex)) and not isinstance(coord, bool)
            # assert isinstance(poly)
        
        return PolygonSpec(
            polygon_vertices=polygon_vertices,
        )
    
    def __init__(self,polygon_vertices):
        
        self._polygon_vertices=polygon_vertices
    
    def polygon_vertices(self):
        return self._polygon_vertices
    
    def to_json(self):
        vertices_copied = [list(v) for v in self._polygon_vertices]
        return {
            "polygon_vertices": vertices_copied,
        }

class EFieldSpecConstantPolygonRegion:
    field_type = "constant_field_polygon_bounded"
    @classmethod
    def from_json(cls, jobj):
        assert jobj["field_type"] == cls.field_type
        
        E_x = jobj["E_x"]
        E_y = jobj["E_y"]
        
        zone_bounds = PolygonSpec.from_json(jobj["zone_bounds"])
        
        return  EFieldSpecConstantPolygonRegion(
            E_x=E_x,
            E_y=E_y,
            zone_bounds=zone_bounds,
        )
    
    def __init__(self, E_x,E_y, zone_bounds):
        """ The E_x,E_y describe the constant field
        
            bounds describe the zone in which the field is on
        """
        self._E_x = E_x
        self._E_y = E_y
        self._zone_bounds = zone_bounds
    
    def zone_bounds(self):
        return self._zone_bounds
    
    def E_x(self):
        return self._E_x
    
    def E_y(self):
        return self._E_y
    
    def to_json(self):
        return {
            "field_type": self.field_type,
            "E_x": self._E_x,
            "E_y": self._E_y,
            "zone_bounds": self._zone_bounds.to_json(),
        }

class ElectricFieldSpec:
    @staticmethod
    def from_json(jobj):
        field_type = jobj["field_type"]
        
        if field_type == EFieldSpecConstantPolygonRegion.field_type:
            return EFieldSpecConstantPolygonRegion.from_json(jobj)
        else:
            raise Exception("Unknown electric field type {}".format(field_type))


class CellForce:
    target_type = "cells"
    @classmethod
    def from_json(cls, jobj):
        assert jobj["target"] == cls.target_type
        if jobj["type"] == CellAreaForce.force_type:
            return CellAreaForce.from_json(jobj)
        elif jobj["type"] == CellPerimeterForce.force_type:
            return CellPerimeterForce.from_json(jobj)
        elif jobj["type"] == UniformElectricForceOnCellBoundary.force_type:
            return UniformElectricForceOnCellBoundary.from_json(jobj)
        elif jobj["type"] == PixelatedElectricForceOnCellBoundary.force_type:
            return PixelatedElectricForceOnCellBoundary.from_json(jobj)
        else:
            raise ValueError("Unknown cell force type '{}'".format(jobj["type"]))

class UniformElectricForceOnCellBoundary(CellForce):
    force_type = "uniform_electric_cell_boundary_force"
    @classmethod
    def from_json(cls, jobj):
        assert jobj["type"] == cls.force_type
        
        return UniformElectricForceOnCellBoundary(
            ElectricFieldSpec.from_json(jobj["field_spec"])
        )
    
    def __init__(self, electric_field_spec):
        self._electric_field_spec = electric_field_spec
    
    def electric_field_spec(self):
        return self._electric_field_spec
    
    def to_json(self):
        return {
            "type": self.force_type,
            "target": self.target_type,
            
            "field_spec": self._electric_field_spec.to_json(),
        }

class PixelatedFieldSpec:
    field_type = "pixelated_electric_field_spec"
    @classmethod
    def from_json(cls, jobj):
        assert jobj["type"] == cls.field_type
        
        grid_origin_x = jobj["grid_origin_x"]
        grid_origin_y = jobj["grid_origin_y"]
        grid_spacing =  jobj["grid_spacing"]
        grid_ncells_x = jobj["grid_ncells_x"]
        grid_ncells_y = jobj["grid_ncells_y"]
        
        assert isinstance(grid_origin_x, float) or isinstance(grid_origin_x, int)
        assert isinstance(grid_origin_y, float) or isinstance(grid_origin_y, int)
        assert isinstance(grid_spacing, float) or isinstance(grid_spacing, int)
        
        assert isinstance(grid_ncells_x, int)
        assert isinstance(grid_ncells_y, int)
        
        field_data = cls.decode_field_data(
            jobj["field_data"],
            grid_ncells_x=grid_ncells_x,
            grid_ncells_y=grid_ncells_y,
        )
        
        return PixelatedFieldSpec(
            grid_origin_x=grid_origin_x,
            grid_origin_y=grid_origin_y,
            grid_spacing=grid_spacing,
            grid_ncells_x=grid_ncells_x,
            grid_ncells_y=grid_ncells_y,
            field_data=field_data,
        )
    
    @classmethod
    def decode_field_data(cls, field_data_arr, grid_ncells_x, grid_ncells_y):
        assert isinstance(field_data_arr, list)
        assert len(field_data_arr) == grid_ncells_x
        for col in field_data_arr:
            assert isinstance(col, list)
            assert len(col) == grid_ncells_y
            
            for field_xy_pair in col:
                assert isinstance(field_xy_pair, list)
                
                assert len(field_xy_pair) == 2 # x and y components of the field at this point
                for el in field_xy_pair:
                    assert isinstance(el, float) or isinstance(el, int) # should be numbers
        
        return copy.deepcopy(field_data_arr)
    
    def __init__(
        self,
        grid_origin_x,
        grid_origin_y,
        grid_spacing,
        grid_ncells_x,
        grid_ncells_y,
        field_data,
        ):
        self._grid_origin_x = grid_origin_x
        self._grid_origin_y = grid_origin_y
        self._grid_spacing = grid_spacing
        self._grid_ncells_x = grid_ncells_x
        self._grid_ncells_y = grid_ncells_y
        self._field_data = field_data
        
    def grid_origin_x(self):
        return self._grid_origin_x
    
    def grid_origin_y(self):
        return self._grid_origin_y
    
    def grid_spacing(self):
        return self._grid_spacing
    
    def grid_ncells_x(self):
        return self._grid_ncells_x
    
    def grid_ncells_y(self):
        return self._grid_ncells_y
    
    def field_data(self):
        return self._field_data
    
    def to_json(self):
        return {
            "type":          self.field_type,
            "grid_origin_x": self._grid_origin_x,
            "grid_origin_y": self._grid_origin_y,
            "grid_spacing":  self._grid_spacing,
            "grid_ncells_x": self._grid_ncells_x,
            "grid_ncells_y": self._grid_ncells_y,
            "field_data": copy.deepcopy(self._field_data),
        }
        
        
class PixElecForceCellParam:
    @classmethod
    def from_json(cls, jobj):
        charge = jobj["charge"]
        assert isinstance(charge, float) or isinstance(charge, int)
        return PixElecForceCellParam(
            charge=charge,
        )
    
    def __init__(self, charge):
        self._charge = charge
        
    def charge(self):
        return self._charge
        
    def to_json(self):
        return {
            "charge": self._charge,
        }

class PixelatedElectricForceOnCellBoundary(CellForce):
    force_type = "pixelated_electric_cell_boundary_force"
    
    @classmethod
    def from_json(cls, jobj):
        assert jobj["type"] == cls.force_type
        
        field_spec = PixelatedFieldSpec.from_json(jobj["field_spec"])
        
        cell_params_jobj = jobj["cell_params"]
        assert isinstance(cell_params_jobj, dict)
        
        
        cell_params_res = {}
        for cell_id, cell_p_jobj in cell_params_jobj.items():
            cell_p_res = PixElecForceCellParam.from_json(cell_p_jobj)
            cell_params_res[cell_id] = cell_p_res
        
        
        return PixelatedElectricForceOnCellBoundary(
            cell_params=cell_params_res,
            field_spec=field_spec,            
        )
    
    def __init__(
        self,
        cell_params,
        field_spec,
    ):
        self._cell_params = cell_params
        self._field_spec = field_spec
    
    def cell_params(self):
        return self._cell_params
    
    def field_spec(self):
        return self._field_spec
    
    def to_json(self):
        cell_params_jobj = {}
        for cell_id, cell_p in self._cell_params.items():
            cell_params_jobj[cell_id] = cell_p.to_json()
        
        
        return {
            "type":        self.force_type,
            "target":      self.target_type,
            "field_spec":  self._field_spec.to_json(),
            "cell_params": cell_params_jobj,
        }

class CellAreaForce(CellForce):
    force_type = "area"
    
    @classmethod
    def from_json(cls, jobj):
        assert jobj["type"] == cls.force_type
        
        return CellAreaForce(
            A0=jobj["A0"],
            kappa=jobj["kappa"],
        )
    
    def __init__(self, A0, kappa):
        self._A0 = A0
        self._kappa = kappa
    
    def A0(self):
        return self._A0
    
    def kappa(self):
        return self._kappa
    
    def to_json(self):
        return {
            "type": self.force_type,
            "target": self.target_type,
            "A0": self._A0,
            "kappa": self._kappa,
        }
    
class CellPerimeterForce(CellForce):
    force_type = "perimeter"
    
    @classmethod
    def from_json(cls, jobj):
        assert jobj["type"] == cls.force_type
        
        return CellPerimeterForce(
            gamma=jobj["gamma"],
            lam=jobj["lambda"],
        )
    
    def __init__(self, gamma, lam):
        self._gamma = gamma
        self._lambda = lam
    
    def gamma(self):
        return self._gamma
    
    def lam(self):
        return self._lambda
    
    def to_json(self):
        return {
            "type": self.force_type,
            "target": self.target_type,
            
            "gamma": self._gamma,
            "lambda": self._lambda,
        }



class VertexForce:
    target_type = "vertices"
    @classmethod
    def from_json(cls, jobj):
        assert jobj["target"] == cls.target_type
        
        if jobj["type"] == ConstantVertexForce.force_type:
            return ConstantVertexForce.from_json(jobj)
        else:
            raise ValueError("Unknown vertex force type '{}'".format(jobj["type"]))

class ConstantVertexForce(VertexForce):
    force_type = "constant_force"
    
    @classmethod
    def from_json(cls, jobj):
        assert jobj["type"] == cls.force_type
        return ConstantVertexForce(
            f_x=jobj["f_x"],
            f_y=jobj["f_y"],
        )
    
    def __init__(self, f_x, f_y):
        self._f_x = f_x
        self._f_y = f_y
    
    def f_x(self):
        return self._f_x
    
    def f_y(self):
        return self._f_y
        
    
    def to_json(self):
        return {
            "type": self.force_type,
            "target": self.target_type,
            
            "f_x": self._f_x,
            "f_y": self._f_y,
        }



class TissueForce:
    @staticmethod
    def from_json(jobj):
        if jobj["target"] == VertexForce.target_type:
            return VertexForce.from_json(jobj)
        elif jobj["target"] == CellForce.target_type:
            return CellForce.from_json(jobj)

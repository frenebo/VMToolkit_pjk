
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
        elif jobj["type"] == ElectricForceOnCellBoundary.force_type:
            return ElectricForceOnCellBoundary.from_json(jobj)
        else:
            raise ValueError("Unknown cell force type '{}'".format(jobj["type"]))

class ElectricForceOnCellBoundary(CellForce):
    force_type = "electric_cell_boundary_force"
    @classmethod
    def from_json(cls, jobj):
        assert jobj["type"] == cls.force_type
        
        return ElectricForceOnCellBoundary(
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



class CellTopology:
    @staticmethod
    def from_json(jobj):
        vertex_ids = jobj["vertex_ids"]
        return CellTopology(
            vertices=vertex_ids,
            is_outer=jobj["is_outer"]
        )
        
    def __init__(self, vertices, is_outer):
        """
        vertices - a list of strings, the ids of vertices belonging to this cell
        is_outer - boolean, whether or not this 'cell' is a boundary cell that exists for sake of 
                   half-edge model
        """
        self._vertex_ids=vertices
        self._is_outer=is_outer
    
    def vertex_ids(self):
        return self._vertex_ids
    
    def is_outer(self):
        return self._is_outer
    
    
    def to_json(self):
        return {
            "vertex_ids": list(self._vertex_ids),
            "is_outer": self._is_outer,
        }

class VertexTopology:
    @staticmethod
    def from_json(jobj):
        return VertexTopology(
            is_boundary=jobj["is_boundary"],
            # num_faces=jobj["num_faces"],
        )
    
    def __init__(self, is_boundary):
        self._is_boundary = is_boundary
        # pass
    
    def is_boundary(self):
        return self._is_boundary
    
    # def num_faces(self):
    #     return self._num_faces
    
    def to_json(self):
        return {
            "is_boundary": self._is_boundary,
        }

class TissueTopology:
    @staticmethod
    def from_json(jobj):
        cell_topologies_parsed = {}
        
        for cell_id, cell_top_json in jobj["cells"].items():
            cell_topologies_parsed[cell_id] = CellTopology.from_json(cell_top_json)
        
        vertex_topologies_parsed = {}
        # if "vertices": not in Tissue
        for vtx_id, vtx_top_json in jobj["vertices"].items():
            vertex_topologies_parsed[vtx_id] = VertexTopology.from_json(vtx_top_json)
        
        # @TODO turn these assertions into checks that give useful exception types and messages that can be handled in a meaningful way
        # assert isinstance(jobj, dict)
        # assert "cells" in jobj
        # assert isinstance(jobj["cells"], dict)
        
        # assert "vertex_ids" in jobj
        # assert isinstance(jobj["vertex_ids"], list)
        
        # for c in jobj["cells"]:
        #     assert isinstance(c, dict)
        #     # assert
        #     assert "verti"
        return TissueTopology(
            cell_topologies=cell_topologies_parsed,
            vertex_topologies=vertex_topologies_parsed
        )
    
    def cell_topologies(self):
        return self._cell_topologies
    
    def vertex_topologies(self):
        return self._vertex_topologies
        
    def __init__(self, cell_topologies, vertex_topologies):
        self._cell_topologies = cell_topologies
        self._vertex_topologies = vertex_topologies
    
    
    def to_json(self):
        json_cell_tops = {}
        for cell_id, cell_top in self._cell_topologies.items():
            json_cell_tops[cell_id] = cell_top.to_json()
        
        json_vtx_tops = {}
        for vtx_id, vtx_top in self._vertex_topologies.items():
            json_vtx_tops[vtx_id] = vtx_top.to_json()
        return {
            "cells": json_cell_tops,
            "vertices": json_vtx_tops,
        }
        

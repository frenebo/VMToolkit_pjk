

class GenericModReq:
    string_properties_names=[]
    str_list_properties_names = []
    typename=None
    @classmethod
    def from_json(cls, jobj):
        processed_args == {}
        for str_prop_name in cls.string_properties_names:
            assert isinstance(jobj[str_prop_name], str)
            processed_args[str_prop_name] = jobj[str_prop_name]
        
        for str_list_pname in cl.str_list_properties_names:
            pval = jobj[str_list_pname]
            assert isinstance(pval, list)
            
            
            for el in pval:
                assert isinstance(el, str)
            processed_args[str_list_pname] = [el for el in pval]
        
        return cls(processed_args)
    
    def __init__(self, properties):
        self._properties = properties
    
    def p(self, prop_name):
        """ Get value of property with the given name
        """
        return self._properties[prop_name]
    
    def to_json(self):
        if self.typename is None:
            raise Exception("Typename must be defined for class")
        jobj = {
            "type": self.typename,
            
        }
        
        for str_prop_name in self.string_properties_names:
            jobj[str_prop_name] = self._properties[str_prop_name]
        for str_list_pname in self.str_list_properties_names:
            jobj[str_list_pname] = [el for el in self._properties[str_list_pname]]
        
        return jobj



class ModReqCreateCellGroup:
    typename = "create_cell_groupf"
    string_properties_names = ["group_name"]

class ModReqCreateVertexGroup(GenericModReq):
    typename = "create_vertex_group"
    string_properties_names = ["group_name"]

class ModReqAddCellsToCellGroup:
    typename = "add_cells_to_cell_group"
    string_properties_names = ["group_name"]
    str_list_properties_names = ["cell_ids"]
    

class ModReqAddVerticesToVertexGroup:
    typename = "add_vertices_to_vertex_group"
    string_properties_names = ["group_name"]
    str_list_properties_names = ["vertex_ids"]

class ModReqRemoveCellsFromCellGroup:
    typename = "remove_cells_from_cell_group"
    string_properties_names = ["group_name"]
    str_list_properties_names = ["cell_ids"]

class ModReqRemoveVerticesFromVertexGroup:
    typename = "remove_vertice_from_vertex_group"
    string_properties_names = ["group_name"]
    str_list_properties_names = ["vertex_ids"]

# class ModReqCreatePerim

# class ModReqCreateForce:
#     typename = "create_force"
    

# class ModReqAddForceToGroup:
#     typename = "add_force_to_group"
#     string_properties_names = ["force_type", "force_id", "target_type", "group_name"]

# class ModReqRemoveForceFromGroup:
#     typename = "remove_force_from_group"
#     # @TODO create force ids
#     string_properties_names = ["force_id", "target_type", "group_name"]

# @TODO add framework for a request parser to branch out in differnt directions...
# class ModReqChangeSimulationIntegratorSetting:
#     typename == "change_simulation_integrator_setting"
#     string_properties_names = ["param_name"]
    

class ModelChangeRequest:
    @staticmethod
    def from_json(jobj):
        all_changereq_classe = [
            ModReqCreateCellGroup,
            ModReqCreateVertexGroup,
            ModReqAddCellsToCellGroup,
            ModReqAddVerticesToVertexGroup,
            ModReqRemoveCellsFromCellGroup,
            ModReqRemoveVerticesFromVertexGroup,
            
            ModReqAddForceToGroup,
            ModReqRemoveForceFromGroup,
            # @TODO add id to forces so that they can be addressed and removed
            
            # ModReq
            ModReqChangeSimulationIntegratorSetting,
            # ModRe
            # ModReqAddForceToVertexGroup,
            # ModReqAddForceToCellGroup,
            
            # ModReqRemoveForceFromVertexGroup,
            # ModReqRemoveForceFromVertexGroup,
            
            # ModChangeFo
        ]
        if jobj["type"] == ModReqCreateCellGroup.typename:
            raise NotImplementedError()
        elif jobj["type"] == ModReqCreateVertexGroup.typename:
            raise NotImplementedError()


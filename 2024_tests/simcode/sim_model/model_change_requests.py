

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
    
    def attr(self, prop_name):
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


class ModReqAddCellsToCellGroup:
    typename = "add_cells_to_cell_group"
    string_properties_names = []
    str_list_properties_names = ["cell_ids"]

class ModReqCreateCellGroup:
    typename = "create_cell_group"
    @staticmethod
    def from_json(jobj):
        return ModReqCreateCellGroup(
            group_name=jobj["group_name"],
        )
    
    def __init__(self, group_name):
        self._group_name = cell_group_name
    
    def group_name(self):
        return self._group_name
    
    def to_json(self):
        return {
            "type": ModReqCreateCellGroup.typename,
            "group_name": self._group_name,
        }

class ModReqCreateVertexGroup(GenericModReq):
    typename == "create_vertex_group"
    string_properties_names = ["group_name"]

class ModelChangeRequest:
    @staticmethod
    def from_json(jobj):
        if jobj["type"] == ModReqCreateCellGroup.typename:
            raise NotImplementedError()
        elif jobj["type"] == ModReqCreateVertexGroup.typename:


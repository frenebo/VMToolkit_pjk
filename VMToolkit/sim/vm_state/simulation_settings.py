


class TopologySettings:
    @staticmethod
    def from_json(jobj):
        t1_sets = T1TransitionSettings.from_json(jobj["T1_transition"])
        
        return TopologySettings(
            T1_transition_settings=t1_sets,
        )
    
    def __init__(self, T1_transition_settings):
        self._T1_transition_settings = T1_transition_settings
    
    def T1_transition_settings(self):
        return self._T1_transition_settings
    
    def to_json(self):
        return {
            "T1_transition": self._T1_transition_settings.to_json(),
        }



class T1TransitionSettings:
    @staticmethod
    def from_json(jobj):
        # if 
        T1_enabled = jobj["enabled"]
        
        if not T1_enabled:
            return T1TransitionSettings(
                enabled=False,
                min_edge_len=None,
                new_edge_len=None,
            )
        else:
            return T1TransitionSettings(
                enabled=True,
                min_edge_len=jobj["min_edge_len"],
                new_edge_len=jobj["new_edge_len"],
            )
        
    def __init__(self, enabled, min_edge_len, new_edge_len):
        self._enabled = enabled
        self._min_edge_len = min_edge_len
        self._new_edge_len = new_edge_len
    
    def enabled(self):
        return self._enabled
    
    def min_edge_len(self):
        return self._min_edge_len
    
    def new_edge_len(self):
        return self._new_edge_len
    
    def to_json(self):
        if not self._enabled:
            return {
                "enabled": False,
                "min_edge_len": self._min_edge_len,
                "new_edge_len": self._new_edge_len,
            }
        else:
            return {
                "enabled": True,
                "min_edge_len": self._min_edge_len,
                "new_edge_len": self._new_edge_len,
            }


class IntegratorSettings:
    @staticmethod
    def from_json(jobj):
        if jobj["integrator_type"] == IntegratorSettingsDormandPrinceRungeKutta.integrator_type:
            return IntegratorSettingsDormandPrinceRungeKutta.from_json(jobj)
        else:
            raise ValueError("Unknown integrator type {}".format(jobj["integrator_type"]))


class IntegratorSettingsDormandPrinceRungeKutta(IntegratorSettings):
    integrator_type = "dormand_prince_runge_kutta"
    
    @classmethod
    def from_json(cls, jobj):
        assert jobj["integrator_type"] == cls.integrator_type
        
        return IntegratorSettingsDormandPrinceRungeKutta(
            vertex_friction_gamma=jobj["vertex_friction_gamma"],
            init_dt=jobj["init_dt"],
            displacement_error_max=jobj["displacement_error_max"],
        )
    
    def __init__(self, vertex_friction_gamma, init_dt, displacement_error_max):
        self._vertex_friction_gamma = vertex_friction_gamma
        self._init_dt = init_dt
        self._displacement_error_max = displacement_error_max
    
    def step_dt(self):
        return self._step_dt
    
    def vertex_friction_gamma(self):
        return self._vertex_friction_gamma
    
    def init_dt(self):
        return self._init_dt
    
    def displacement_error_max(self):
        return self._displacement_error_max
    
    def to_json(self):
        return {
            "integrator_type": self.integrator_type,
            "init_dt": self._init_dt,
            "displacement_error_max": self._displacement_error_max,
            "vertex_friction_gamma": self._vertex_friction_gamma
        }



class SimulationSettings:
    @staticmethod
    def from_json(jobj):
        return SimulationSettings(
            integrator_settings=IntegratorSettings.from_json(jobj["integrator"]),
            topology_settings=TopologySettings.from_json(jobj["topology_settings"]),
        )
    
    def __init__(self, integrator_settings, topology_settings):
        self._integrator_settings = integrator_settings
        self._topology_settings = topology_settings
        
    def integrator_settings(self):
        return self._integrator_settings
    
    def topology_settings(self):
        return self._topology_settings
    
    def to_json(self):
        return {
            "integrator": self._integrator_settings.to_json(),
            "topology_settings": self._topology_settings.to_json(),
        }
    
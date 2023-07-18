import numpy as np
from idmtools.entities.simulation import Simulation
from typing import Dict, Any


##################################################
# Sweeping utility functions
##################################################
def set_param(simulation: Simulation, param: str, value: Any) -> Dict[str, Any]:
    """
    Set specific parameter value
    Args:
        simulation: idmtools Simulation
        param: parameter
        value: new value

    Returns:
        dict
    """
    return simulation.task.set_parameter(param, value)


class ItvFn:
    """
    Requirements:
     - func is a method that takes campaign as first parameter
     - func return a dict
     - if original is just returning an event, we need to make a wrapper
    """

    def __init__(self, func, *args, **kwargs):
        self.func = func
        self.args = args
        self.kwargs = kwargs

    def __call__(self, simulation: Simulation):
        import emod_api.campaign as campaign
        campaign.reset()

        md = self.func(campaign, *self.args, **self.kwargs)

        # add events
        events = campaign.campaign_dict["Events"]
        simulation.task.campaign.add_events(events)

        # update config for adhoc events
        adhoc_events = campaign.get_adhocs()
        if len(adhoc_events) > 0:
            print("Found adhoc events in campaign. Needs some special processing behind the scenes.")
            if "Custom_Individual_Events" in simulation.task.config.parameters:
                ev_exist = set(simulation.task.config.parameters.Custom_Individual_Events)
                ev_addhoc = set(adhoc_events.keys())
                simulation.task.config.parameters.Custom_Individual_Events = list(ev_exist.union(ev_addhoc))
            else:
                simulation.task.config.parameters.Custom_Individual_Events = list(set(adhoc_events.keys()))

        # Make sure we cast numpy types into normal system types
        if md:
            for k, v in md.items():
                if isinstance(v, (np.int64, np.float64, np.float32, np.uint32, np.int16, np.int32)):
                    md[k] = v.item()

        return md


class CfgFn:
    """
    Requirements:
     - func is a method that takes config as first parameter
     - func return a dict
    """

    def __init__(self, func, *args, **kwargs):
        self.func = func
        self.args = args
        self.kwargs = kwargs

    def __call__(self, simulation: Simulation):
        md = self.func(simulation.task.config, *self.args, **self.kwargs)

        # Make sure we cast numpy types into normal system types
        if md:
            for k, v in md.items():
                if isinstance(v, (np.int64, np.float64, np.float32, np.uint32, np.int16, np.int32)):
                    md[k] = v.item()

        return md


class SwpFn:
    """
    This works for sweeping on report, demographics, migrations and climate
    Requirements:
     - func is a method that takes task as first parameter
     - func return a dict
    """

    def __init__(self, func, *args, **kwargs):
        self.func = func
        self.args = args
        self.kwargs = kwargs

    def __call__(self, simulation: Simulation):
        md = self.func(simulation.task, *self.args, **self.kwargs)

        # Make sure we cast numpy types into normal system types
        if md:
            for k, v in md.items():
                if isinstance(v, (np.int64, np.float64, np.float32, np.uint32, np.int16, np.int32)):
                    md[k] = v.item()

        return md

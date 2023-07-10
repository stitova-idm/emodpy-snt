from typing import Any, Dict
from functools import partial
from idmtools.builders import SimulationBuilder
from idmtools.entities.simulation import Simulation

import manifest  # required in helper functions
import params

platform = None


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


###################################
# Common interface
###################################
def get_sweep_builders(**kwargs):
    global platform
    platform = kwargs.get('platform', None)

    builder = SimulationBuilder()

    # Test
    if params.test_run:
        builder.add_sweep_definition(partial(set_param, param='Run_Number'), range(params.num_seeds))
        print(builder.count)
    else:
        pass

    return [builder]

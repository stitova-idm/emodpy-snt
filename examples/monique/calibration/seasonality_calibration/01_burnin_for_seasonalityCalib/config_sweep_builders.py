import params
from functools import partial
from idmtools.builders import SimulationBuilder
from snt.utility.sweeping import set_param, ItvFn, CfgFn, sweep_functions

platform = None


###################################
# Common interface
###################################
def get_sweep_builders(**kwargs):
    """
    Build simulation builders.
    Args:
        kwargs: User inputs may overwrite the entries in the block.

    Returns:
        lis of Simulation builders
    """
    global platform
    platform = kwargs.get('platform', None)
    builder = SimulationBuilder()

    funcs_list = [[partial(set_param, param='Run_Number', value=x),
                   partial(set_param, param='MaxHab', value=params.max_habitat_value),
                   partial(set_param, param='ConstHab', value=params.const_habitat),
                   partial(set_param, param='Admin_Name', value=my_rep_admin),
                   ]
                  for my_rep_admin in [params.rep_admin]
                  for x in range(params.num_seeds)
                  ]

    builder.add_sweep_definition(sweep_functions, funcs_list)

    return [builder]

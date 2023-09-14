import os
import pandas as pd
from tqdm import tqdm
from typing import List
from functools import partial
from idmtools.builders import SimulationBuilder
from idmtools.entities.simulation import Simulation
from emodpy_malaria.reporters.builtin import add_report_event_counter
from snt.utility.sweeping import set_param, ItvFn, CfgFn
from snt.hbhi.set_up_general import setup_ds

import manifest
import params

platform = None


def sweep_interventions(simulation: Simulation, func_list: List):
    # Specially handel add_all_interventions to add report!
    tags_updated = {}
    for func in func_list:
        tags = func(simulation)
        if tags:
            if isinstance(func, ItvFn):
                fname = func.func.__name__
                if fname == 'add_all_interventions':
                    add_report_event_counter(simulation.task, manifest, event_trigger_list=tags["events"])
            else:
                tags_updated.update(tags)
    return tags_updated


###################################
# Common interface
###################################
def get_sweep_builders(**kwargs):
    global platform
    platform = kwargs.get('platform', None)

    builder = SimulationBuilder()

    # Important DFs
    master_df = params.master_df
    lhdf = pd.read_csv(os.path.join(params.iopath, params.larval_hab_csv))
    rel_abund_df = pd.read_csv(os.path.join(params.iopath, params.rel_abund_csv))
    rel_abund_df = rel_abund_df.set_index('DS_Name')

    samp_df = params.samp_df
    ds_list = params.ds_list

    # Test
    if params.test_run:
        builder.add_sweep_definition(partial(set_param, param='Run_Number'), range(params.num_seeds))
        return [builder]

    # BUILDER
    int_sweeps = []
    for my_ds in tqdm(ds_list):
        samp_ds = samp_df[samp_df.archetype == my_ds]
        # For each archetype, each Habitat_Multiplier has same reduce_id and seed1 by design of preprocessing
        samp_ds1 = samp_ds[['Habitat_Multiplier', 'reduce_id', 'seed1']].drop_duplicates()
        for r, row in samp_ds1.iterrows():
            int_funcs = []
            cnf = CfgFn(setup_ds,
                        manifest=manifest,
                        platform=platform,
                        my_ds=my_ds,
                        archetype_ds=master_df.at[my_ds, 'archetype'],
                        pull_from_serialization=params.pull_from_serialization,
                        burnin_id='',
                        ser_date=0,
                        rel_abund_df=rel_abund_df,
                        lhdf=lhdf,
                        demographic_suffix=params.demographic_suffix,
                        climate_prefix=params.climate_prefix,
                        climate_suffix=params.climate_suffix,
                        use_arch_burnin=params.use_arch_burnin,
                        use_arch_input=params.use_arch_input,
                        hab_multiplier=row['Habitat_Multiplier'])
            int_funcs.append(cnf)
            int_funcs.append(partial(set_param, param='Habitat_Multiplier',
                                     value=row['Habitat_Multiplier']))
            int_funcs.append(partial(set_param, param='Run_Number', value=row['seed1']))
            int_funcs.append(
                partial(set_param, param='Sample_ID', value=row['reduce_id']))

            int_sweeps.append(int_funcs)

    builder.add_sweep_definition(sweep_interventions, int_sweeps)
    print(builder.count)

    return [builder]

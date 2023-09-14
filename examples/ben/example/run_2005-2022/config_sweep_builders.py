import os
import pandas as pd
from tqdm import tqdm
from typing import List
from functools import partial
from scipy import interpolate
from idmtools.builders import SimulationBuilder
from idmtools.entities.simulation import Simulation
from emodpy_malaria.reporters.builtin import add_report_event_counter
from snt.utility.sweeping import set_param, ItvFn, CfgFn, SwpFn
from snt.hbhi.set_up_general import setup_ds
from snt.hbhi.set_up_interventions import add_all_interventions, update_smc_access_ips
from snt.hbhi.utils import tryread_df

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
                    if len(tags["events"]) > 0 and params.add_event_report:
                        add_report_event_counter(simulation.task, manifest, event_trigger_list=tags["events"],
                                                 start_day=0, end_day=params.years * 365)
            else:
                tags_updated.update(tags)
    return tags_updated


# BUILDER
def lin_interpolate(years, vals):
    val_interp = interpolate.interp1d(years, vals)
    full_vals = val_interp(range(years[0], years[-1] + 1))
    return (full_vals)


###################################
# Common interface
###################################
def get_sweep_builders(**kwargs):
    global platform
    platform = kwargs.get('platform', None)

    builder = SimulationBuilder()

    # Treatment-seeking
    hs_df = tryread_df(os.path.join(params.iopath, 'simulation_inputs',
                                    'scenarios', 'cm_2005-2022.csv'))
    # ITNs
    itn_df = tryread_df(os.path.join(params.iopath, 'simulation_inputs',
                                     'scenarios', 'itn_2005-2022.csv'))
    # SMC
    smc_df = tryread_df(os.path.join(params.iopath, 'simulation_inputs',
                                     'scenarios', 'smc_2005-2022.csv'))
    smc_df['high_access_5_10'] = smc_df['high_access_O5']  # update_smc_access_ips hard-coded the column name...
    smc_df['coverage_low_access_U5'] = 0
    smc_df['coverage_low_access_O5'] = 0

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

    int_suite = params.int_suite

    # BUILDER
    int_sweeps = []
    for my_ds in tqdm(ds_list):
        samp_ds = samp_df[samp_df.DS_Name == my_ds].copy()
        for r, row in samp_ds.iterrows():
            hs_ds = hs_df.copy()[hs_df[int_suite.hs_ds_col] == my_ds]
            full_covs = lin_interpolate([2005, 2013, 2017, 2020, 2022],
                                        [0, row['CM_cov_2013'], row['CM_cov_2017'],
                                         row['CM_cov_2020'], row['CM_cov_2020']])
            hs_ds['U5_coverage'] = full_covs * hs_ds['U5_coverage']
            hs_ds['adult_coverage'] = full_covs * hs_ds['adult_coverage']

            itn_ds = itn_df.copy()[itn_df[int_suite.itn_ds_col] == my_ds]
            itn_ds = itn_ds[itn_ds['year'] <= 2022]
            itn_ds = itn_ds.reset_index()

            # TODO: Make code below more defensive
            itn_ds.loc[0, int_suite.itn_cov_cols] = itn_ds.loc[0, int_suite.itn_cov_cols] * row['ITN_2011']
            itn_ds.loc[1, int_suite.itn_cov_cols] = itn_ds.loc[1, int_suite.itn_cov_cols] * row['ITN_2014']
            itn_ds.loc[2, int_suite.itn_cov_cols] = itn_ds.loc[2, int_suite.itn_cov_cols] * row['ITN_2017']
            itn_ds.loc[3, int_suite.itn_cov_cols] = itn_ds.loc[3, int_suite.itn_cov_cols] * row['ITN_2020']

            cnf = CfgFn(setup_ds,
                        manifest=manifest,
                        platform=platform,
                        my_ds=my_ds,
                        archetype_ds=master_df.at[my_ds, 'archetype'],
                        pull_from_serialization=params.pull_from_serialization,
                        burnin_id=params.burnin_id,
                        ser_date=params.ser_date,
                        rel_abund_df=rel_abund_df,
                        lhdf=lhdf,
                        demographic_suffix=params.demographic_suffix,
                        climate_prefix=params.climate_prefix,
                        climate_suffix=params.climate_suffix,
                        use_arch_burnin=params.use_arch_burnin,
                        use_arch_input=params.use_arch_input,
                        hab_multiplier=row['Habitat_Multiplier'],
                        serialize_match_tag=['Habitat_Multiplier'],
                        serialize_match_val=[float(row['Habitat_Multiplier'])])
            int_f = ItvFn(add_all_interventions,
                          int_suite=int_suite,
                          my_ds=my_ds,
                          hs_df=hs_ds,
                          itn_df=itn_ds,
                          smc_df=smc_df,
                          addtl_smc_func=update_smc_access_ips  # Change IP every year
                          )

            for x in range(params.num_seeds):
                int_funcs = []

                int_funcs.append(cnf)
                int_funcs.append(int_f)
                int_funcs.append(partial(set_param, param='DS_Name', value=my_ds))
                int_funcs.append(partial(set_param, param='Sample_ID', value=row['id']))
                int_funcs.append(partial(set_param, param='Run_Number', value=row['seed2'] + x))

                int_sweeps.append(int_funcs)

    builder.add_sweep_definition(sweep_interventions, int_sweeps)
    print(builder.count)

    return [builder]

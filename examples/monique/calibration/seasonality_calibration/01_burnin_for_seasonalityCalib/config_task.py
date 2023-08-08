import os
import pandas as pd
from pathlib import Path
from emodpy import emod_task
from emodpy_malaria.interventions.outbreak import add_outbreak_individual
from snt.helpers_add_interventions import add_all_interventions

import manifest
import params

platform = None


#####################################
# Utility functions
#####################################

def _config_reports(task):
    """
    Add reports

    Args:
        task: EMODTask

    Returns: None

    """
    from emodpy_malaria.reporters.builtin import add_report_malaria_filtered

    add_report_malaria_filtered(task, manifest,
                                start_day=(params.years - 5) * 365,
                                end_day=params.years * 365)


#####################################
# Create EMODTask
#####################################

def build_campaign():
    """
    Adding required interventions

    Returns:
        campaign object
    """

    import emod_api.campaign as campaign
    # passing in schema file to verify that everything is correct.
    campaign.schema_path = manifest.schema_file

    add_outbreak_individual(campaign, start_day=182, demographic_coverage=0.01, timesteps_between_repetitions=365)

    project_path = manifest.project_path
    scen_df = params.scen_df
    scen_index = params.scen_index

    # INTERVENTIONS
    # health-seeking
    if (not pd.isna(scen_df.at[scen_index, 'CM_filename'])) and (not (scen_df.at[scen_index, 'CM_filename'] == '')):
        hs_df = pd.read_csv(
            os.path.join(project_path, 'simulation_inputs', '%s.csv' % scen_df.at[scen_index, 'CM_filename']))
    else:
        hs_df = pd.DataFrame()
    # NMFs
    if (not pd.isna(scen_df.at[scen_index, 'NMF_filename'])) and (not (scen_df.at[scen_index, 'NMF_filename'] == '')):
        nmf_df = pd.read_csv(
            os.path.join(project_path, 'simulation_inputs', '%s.csv' % scen_df.at[scen_index, 'NMF_filename']))
    else:
        nmf_df = pd.DataFrame()
    # ITNs
    if (not pd.isna(scen_df.at[scen_index, 'ITN_filename'])) and (not (scen_df.at[scen_index, 'ITN_filename'] == '')):
        itn_df = pd.read_csv(
            os.path.join(project_path, 'simulation_inputs', '%s.csv' % scen_df.at[scen_index, 'ITN_filename']))
    else:
        itn_df = pd.DataFrame()
    if (not pd.isna(scen_df.at[scen_index, 'ANC_ITN_filename'])) and (
            not (scen_df.at[scen_index, 'ANC_ITN_filename'] == '')):
        itn_anc_df = pd.read_csv(
            os.path.join(project_path, 'simulation_inputs', '%s.csv' % scen_df.at[scen_index, 'ANC_ITN_filename']))
    else:
        itn_anc_df = pd.DataFrame()
    if (not pd.isna(scen_df.at[scen_index, 'EPI_ITN_filename'])) and (
            not (scen_df.at[scen_index, 'EPI_ITN_filename'] == '')):
        itn_epi_df = pd.read_csv(
            os.path.join(project_path, 'simulation_inputs', '%s.csv' % scen_df.at[scen_index, 'EPI_ITN_filename']))
    else:
        itn_epi_df = pd.DataFrame()

    itn_decay_params = pd.read_csv(os.path.join(project_path, 'simulation_inputs', 'itn_discard_decay_params.csv'))
    itn_use_seasonality = pd.read_csv(os.path.join(project_path, 'simulation_inputs', 'ITN_use_seasonality.csv'))

    # IRS
    if (not pd.isna(scen_df.at[scen_index, 'IRS_filename'])) and (not (scen_df.at[scen_index, 'IRS_filename'] == '')):
        irs_df = pd.read_csv(
            os.path.join(project_path, 'simulation_inputs', '%s.csv' % scen_df.at[scen_index, 'IRS_filename']))
    else:
        irs_df = pd.DataFrame()
    # SMC
    if (not pd.isna(scen_df.at[scen_index, 'SMC_filename'])) and (not (scen_df.at[scen_index, 'SMC_filename'] == '')):
        smc_df = pd.read_csv(
            os.path.join(project_path, 'simulation_inputs', '%s.csv' % scen_df.at[scen_index, 'SMC_filename']))
    else:
        smc_df = pd.DataFrame()

    add_all_interventions(campaign,
                          hfca=params.rep_admin,
                          itn_anc_adult_birthday_years=params.itn_anc_adult_birthday_years,
                          itn_decay_params=itn_decay_params,
                          itn_use_seasonality=itn_use_seasonality,
                          hs_df=hs_df,
                          nmf_df=nmf_df,
                          itn_df=itn_df,
                          itn_anc_df=itn_anc_df,
                          itn_epi_df=itn_epi_df,
                          irs_df=irs_df,
                          smc_df=smc_df)

    return campaign


def set_param_fn(config):
    """
    This function is a callback that is passed to emod-api.config to set parameters The Right Way.

    Args:
        config:

    Returns:
        configuration settings
    """

    # You have to set simulation type explicitly before you set other parameters for the simulation
    # sets "default" malaria parameters as determined by the malaria team
    import emodpy_malaria.malaria_config as malaria_config
    config = malaria_config.set_team_defaults(config, manifest)

    from set_config import set_config
    set_config(config)

    return config


def get_task(**kwargs):
    """
    This function is designed to create and config a Task

    Args:
        **kwargs: optional parameters

    Returns:
        task

    """
    global platform
    platform = kwargs.get('platform', None)

    # Create EMODTask
    print("Creating EMODTask...")
    task = emod_task.EMODTask.from_default2(
        config_path=None,
        eradication_path=manifest.eradication_path,
        schema_path=manifest.schema_file,
        campaign_builder=build_campaign,
        param_custom_cb=set_param_fn,
        ep4_custom_cb=None,
    )

    # Config demographics
    demog_path = os.path.join(manifest.input_dir, params.demographics_file)
    task.common_assets.add_asset(demog_path, relative_path=str(Path(params.demographics_file).parent))

    # More stuff to add task, like reports...
    _config_reports(task)

    return task

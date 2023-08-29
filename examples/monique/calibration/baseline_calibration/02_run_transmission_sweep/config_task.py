import os
from pathlib import Path
from emodpy import emod_task
import manifest
import params

platform = None


#####################################
# Utility functions
#####################################

def _config_reports(task):
    """
    Add reports.

    Args:
        task: EMODTask
    Returns:
        None
    """
    from emodpy_malaria.reporters.builtin import add_report_malaria_filtered
    from emodpy_malaria.reporters.builtin import add_malaria_summary_report
    add_report_malaria_filtered(task, manifest,
                                start_day=0,
                                end_day=params.years * 365)

    for year in range(params.years):
        add_malaria_summary_report(task, manifest,
                                   start_day=365 * year,
                                   end_day=365 * year + 365,
                                   age_bins=[0.25, 5, 15, 30, 50, 125],
                                   reporting_interval=30,
                                   parasitemia_bins=[10, 50, 1e9],
                                   filename_suffix='Monthly%d' % (year + params.start_year)
                                   )

    #####################################
    # Create EMODTask
    #####################################


def build_campaign():
    """
    Adding required interventions.

    Returns:
        campaign object
    """

    import emod_api.campaign as campaign
    # passing in schema file to verify that everything is correct.
    campaign.schema_path = manifest.schema_file

    from emodpy_malaria.interventions.outbreak import add_outbreak_individual
    add_outbreak_individual(campaign, start_day=35, demographic_coverage=0.002, timesteps_between_repetitions=73)

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
    This function is designed to create and config a Task.

    Args:
        kwargs: optional parameters
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

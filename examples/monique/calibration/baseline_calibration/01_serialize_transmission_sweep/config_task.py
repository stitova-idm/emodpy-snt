from emodpy import emod_task
from emodpy_malaria.interventions.outbreak import add_outbreak_individual

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
    from emodpy_malaria.reporters.builtin import add_malaria_summary_report
    add_report_malaria_filtered(task, manifest,
                                start_day=(params.years - 3) * 365,
                                end_day=params.years * 365
                                )
    add_malaria_summary_report(task, manifest,
                               start_day=(params.years - 1) * 365,
                               reporting_interval=30,
                               age_bins=[0.25, 5, 15, 125],
                               parasitemia_bins=[10, 50, 1e9]
                               )


def _config_campaign(campaign):
    """
    Add interventions

    Args:
        campaign: emod_api.campaign

    Returns: None

    """
    add_outbreak_individual(campaign, start_day=35, demographic_coverage=0.002, timesteps_between_repetitions=73)
    pass


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

    _config_campaign(campaign)

    return campaign


def set_config_parameters(config):
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
        param_custom_cb=set_config_parameters,
        ep4_custom_cb=None,
    )

    # More stuff to add to task, like reports...
    _config_reports(task)

    return task

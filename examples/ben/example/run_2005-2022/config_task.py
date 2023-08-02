import manifest
import params
from emodpy import emod_task
from emodpy_malaria.interventions.outbreak import add_outbreak_individual

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
    from snt.hbhi.set_up_general import initialize_reports
    initialize_reports(task, manifest,
                       params.event_reporter,
                       params.filtered_report,
                       params.years,
                       params.yr_plusone)


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

    add_outbreak_individual(campaign, start_day=182, demographic_coverage=0.01, repetitions=-1,
                            timesteps_between_repetitions=365)

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

    # Add assets corresponding to the filename parameters set in set_input_files.
    from snt.hbhi.set_up_general import add_input_files
    master_df = params.master_df
    ds_list = params.ds_list

    for my_ds in ds_list:
        add_input_files(task,
                        iopath=params.iopath,
                        my_ds=my_ds,
                        archetype_ds=master_df.at[my_ds, 'archetype'],
                        demographic_suffix=params.demographic_suffix,
                        climate_suffix=params.climate_suffix,
                        climate_prefix=params.climate_prefix,
                        use_archetype=params.use_arch_input
                        )

    # More stuff to add task, like reports...
    _config_reports(task)

    return task

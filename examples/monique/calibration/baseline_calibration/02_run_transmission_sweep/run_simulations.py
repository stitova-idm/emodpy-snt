import manifest
import params
from idmtools.core.platform_factory import Platform
from idmtools.entities.experiment import Experiment
from idmtools.entities.templated_simulation import TemplatedSimulations
from snt.utility.emod_api_utils import suppress_warnings


def _print_params():
    """
    Just a useful convenient function for the user.
    """
    print("expname: ", params.expname)
    print("pull_from_serialization: ", params.pull_from_serialization)
    print("population_size: ", params.population_size)
    print("years: ", params.years)
    print("num_seeds: ", params.num_seeds)
    print("burnin_id: ", params.burnin_id)


def _post_run(experiment: Experiment, **kwargs):
    """
    Add extra work after run experiment.
    Args:
        experiment: idmtools Experiment
        kwargs: additional parameters
    Return:
        None
    """
    # Save experiment id to file to be used by snakefile
    with open(r"monique\\calibration\\baseline_calibration\\02_run_transmission_sweep\\experiment_id.txt", "w") as fd:
        fd.write(experiment.uid.hex)
    pass


def _config_experiment(**kwargs):
    """
    Build experiment from task and builder. task is EMODTask. builder is SimulationBuilder used for config parameter sweeping.
    Args:
        kwargs: additional parameters
    Return:
        experiment
    """
    from config_sweep_builders import get_sweep_builders
    from config_task import get_task

    builders = get_sweep_builders(**kwargs)

    task = get_task(**kwargs)

    if manifest.sif_id:
        task.set_sif(manifest.sif_id)

    ts = TemplatedSimulations(base_task=task, builders=builders)

    experiment = Experiment.from_template(ts, name=params.expname)

    return experiment


def run_experiment(show_warnings: bool = True, **kwargs):
    """
    Get configured calibration and run.
    Args:
        show_warnings: True/False
        kwargs: user inputs
    Returns:
        None
    """
    # make sure pass platform through
    kwargs['platform'] = platform

    # Suppress emod_api warnings
    suppress_warnings(show_warnings=show_warnings)

    _print_params()

    experiment = _config_experiment(**kwargs)
    experiment.run(wait_until_done=True, wait_on_done=False)
    _post_run(experiment, **kwargs)


if __name__ == "__main__":
    platform = Platform('CALCULON', node_group='idm_48cores')
    # platform = Platform('IDMCLOUD', node_group='emod_abcd')

    # If you don't have Eradication, un-comment out the following to download Eradication
    # import emod_malaria.bootstrap as dtk
    #
    # dtk.setup(pathlib.Path(manifest.eradication_path).parent)
    # os.chdir(os.path.dirname(__file__))
    # print("...done.")
    run_experiment(show_warnings=False)

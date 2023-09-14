import warnings
from idmtools.entities.simulation import Simulation
from idmtools.registry.hook_specs import IDMTOOLS_HOOKS
from idmtools.registry.functions import FunctionPluginManager
import emod_api.schema_to_class as s2c
from pluggy import HookimplMarker
from threading import Lock

function_hook_impl = HookimplMarker(IDMTOOLS_HOOKS)

MSG = 'We are hiding warnings just for clarity purpose. Users should pay attention to the warnings and adjust coding to avoid such warnings in the future.'

show_warnings_lock = Lock()
SHOW_WARNINGS = None
WARNINGS_ONCE = None


def initialize_plugins(**kwargs):
    """
    Setup plugins.
    Args:
        kwargs: user inputs
    Returns:
        None
    """
    global WARNINGS_ONCE
    WARNINGS_ONCE = kwargs.get('show_warnings_once', None)
    _initialize_warnings(WARNINGS_ONCE)

    # register plugins
    pm = FunctionPluginManager.instance()
    pm.register(Plugin_pre_create())


def _initialize_warnings(show_warnings_once: bool = None):
    """
    Turn on/off emod_api warnings.
    Args:
        show_warnings_once: None/True/False
    Returns:
        None
    """
    global SHOW_WARNINGS
    if show_warnings_once is None:
        SHOW_WARNINGS = False
        display_info()
    else:
        SHOW_WARNINGS = True
        if show_warnings_once:
            display_info()

    s2c.show_warnings = SHOW_WARNINGS


def _suppress_api_warnings():
    """
    Show warnings only for one simulation.
    Args:
        platform: idmtools platform
    Returns:
        None
    """
    global SHOW_WARNINGS
    if SHOW_WARNINGS and WARNINGS_ONCE:
        show_warnings_lock.acquire()
        SHOW_WARNINGS = False
        show_warnings_lock.release()
    else:
        if s2c.show_warnings and WARNINGS_ONCE:
            show_warnings_lock.acquire()
            s2c.show_warnings = False
            show_warnings_lock.release()


def display_info():
    warnings.warn(MSG, UserWarning)


class Plugin_pre_create:
    """A 2nd hook implementation namespace."""

    @function_hook_impl
    def idmtools_platform_pre_create_item(self, item):
        """
        This callback is called by the pre_create of each object type on a platform.
        An item can be a suite, workitem, simulation, asset collection or an experiment.
        Args:
            item: idmtools entity
        Returns:
            None
        """
        if isinstance(item, Simulation):
            _suppress_api_warnings()

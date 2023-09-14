import time
import functools
from COMPS.Data import Experiment
from COMPS.Data import QueryCriteria


def retry_function(func, wait=1.5, max_retries=5):
    """
    Decorator allowing to retry the call to a function with some time in between.
    Args:
        func:
        wait:
        max_retries:

    Returns:
        None

    Usage:

    @retry_function
    def my_func():
        pass

    @retry_function(max_retries=10, wait=2)
    def my_func():
        pass

    """

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        retExc = None
        for i in range(max_retries):
            try:
                return func(*args, **kwargs)
            except RuntimeError as r:
                # Immediately raise if this is an error.
                # COMPS is reachable so let's be clever and trust COMPS
                if str(r) == "404 NotFound - Failed to retrieve experiment for given id":
                    raise r
            except Exception as e:
                retExc = e
                time.sleep(wait)
        raise retExc if retExc else Exception()

    return wrapper


@retry_function
def get_experiments_by_name(name: str, user=None):
    filters = ["name~{}".format(name)]
    if user:
        filters.append("owner={}".format(user))
    return Experiment.get(query_criteria=QueryCriteria().where(filters))


@retry_function
def get_experiments_by_name(name: str, user=None):
    filters = ["name~{}".format(name)]
    if user:
        filters.append("owner={}".format(user))
    return Experiment.get(query_criteria=QueryCriteria().where(filters))


@retry_function
def get_most_recent_experiment_id_by_name(name: str, user=None):
    experiments = {e.date_created: e for e in get_experiments_by_name(name, user)}
    return experiments[max(experiments.keys())].id

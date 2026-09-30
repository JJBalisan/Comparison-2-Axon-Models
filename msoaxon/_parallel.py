"""Running independent simulations across worker processes."""

import multiprocessing
import os
from concurrent.futures import ProcessPoolExecutor


def pool_context(method=None):
    """The multiprocessing context for worker pools: the platform default, except
    that fork is replaced by forkserver.

    Forking a process that runs threads (a pytest-xdist worker, a numerical library's
    thread pool) can deadlock the child, and Python 3.12+ warns when it happens;
    Python 3.14 makes forkserver the Linux default for this reason. macOS and
    Windows already default to spawn. method: the default to assume (for tests).
    """
    method = method or multiprocessing.get_context().get_start_method()
    return multiprocessing.get_context("forkserver" if method == "fork" else method)


def process_pool(max_workers=None):
    """A ProcessPoolExecutor using pool_context(). Use this for every worker pool."""
    return ProcessPoolExecutor(max_workers=max_workers, mp_context=pool_context())


def map_tasks(fn, tasks, workers=None, executor=None, chunksize=1):
    """[fn(t) for t in tasks], in parallel unless workers == 1.

    executor: an existing pool to run on. Scripts that make several parallel calls
    should open one process_pool() and pass it to each: a fresh pool per call
    re-imports msoaxon in every worker and starts with an empty pre-stimulus cache.
    Otherwise a pool of `workers` processes (default: all CPUs) is made for this call.
    fn must be a module-level function so the workers can import it.
    """
    if executor is not None:
        return list(executor.map(fn, tasks, chunksize=chunksize))
    workers = min(workers or os.cpu_count() or 1, len(tasks))
    if workers <= 1:  # also no tasks at all, which a pool of 0 workers would reject
        return [fn(t) for t in tasks]
    with process_pool(workers) as pool:
        return list(pool.map(fn, tasks, chunksize=chunksize))

"""Running independent simulations across worker processes."""

import multiprocessing
import os
from concurrent import futures

# every module whose functions run in pool workers
_PRELOAD = [
    "msoaxon.coincidence",
    "msoaxon.somatic",
    "msoaxon.threshold",
    "msoaxon.measure",
]


def pool_context(method=None):
    """Returns the multiprocessing context for worker pools.

    That is the platform default, except that fork is replaced by
    forkserver. Forking a process that runs threads (a pytest-xdist
    worker, a numerical library's thread pool) can deadlock the child,
    and Python 3.12+ warns when it happens; Python 3.14 makes forkserver
    the Linux default for this reason. macOS and Windows already default
    to spawn.

    Args:
        method: The default start method to assume (for tests); None
            reads the platform's.

    Returns:
        A multiprocessing context.
    """
    method = method or multiprocessing.get_context().get_start_method()
    if method != "fork":
        return multiprocessing.get_context(method)
    ctx = multiprocessing.get_context("forkserver")
    # the server imports these once, so each worker forks with msoaxon, numpy and
    # scipy already loaded instead of importing them itself (nearly fork's start-up
    # cost; the server runs no other threads, so forking it is safe). No effect once
    # the server is running.
    ctx.set_forkserver_preload(_PRELOAD)
    return ctx


def process_pool(max_workers=None):
    """Returns a ProcessPoolExecutor using pool_context().

    Use this for every worker pool.

    Args:
        max_workers: The number of processes (default: all CPUs).

    Returns:
        A concurrent.futures.ProcessPoolExecutor.
    """
    return futures.ProcessPoolExecutor(
        max_workers=max_workers, mp_context=pool_context()
    )


def map_tasks(fn, tasks, workers=None, executor=None, chunksize=1):
    """Computes [fn(t) for t in tasks], in parallel unless workers == 1.

    Args:
        fn: The function to apply; it must be a module-level function
            so the workers can import it.
        tasks: The arguments, one per call.
        workers: Without an executor, the size of a pool made for this
            call (default: all CPUs); 1 (or a single task) runs
            in-process.
        executor: An existing pool to run on. Scripts that make several
            parallel calls should open one process_pool() and pass it to
            each: a fresh pool per call re-imports msoaxon in every
            worker and starts with an empty pre-stimulus cache.
        chunksize: Tasks sent to a worker at a time.

    Returns:
        The results, in the order of `tasks`.
    """
    if executor is not None:
        return list(executor.map(fn, tasks, chunksize=chunksize))
    workers = min(workers or os.cpu_count() or 1, len(tasks))
    if (
        workers <= 1
    ):  # also no tasks at all, which a pool of 0 workers would reject
        return [fn(t) for t in tasks]
    with process_pool(workers) as pool:
        return list(pool.map(fn, tasks, chunksize=chunksize))

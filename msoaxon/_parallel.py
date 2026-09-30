"""Running independent simulations across worker processes."""

import os
from concurrent.futures import ProcessPoolExecutor


def map_tasks(fn, tasks, workers=None, executor=None, chunksize=1):
    """[fn(t) for t in tasks], in parallel unless workers == 1.

    executor: an existing pool to run on. Scripts that make several parallel calls
    should open one `ProcessPoolExecutor` and pass it to each: a fresh pool per call
    re-imports msoaxon in every worker and starts with an empty pre-stimulus cache.
    Otherwise a pool of `workers` processes (default: all CPUs) is made for this call.
    fn must be a module-level function so the workers can import it.
    """
    if executor is not None:
        return list(executor.map(fn, tasks, chunksize=chunksize))
    workers = min(workers or os.cpu_count() or 1, len(tasks))
    if workers == 1:
        return [fn(t) for t in tasks]
    with ProcessPoolExecutor(max_workers=workers) as pool:
        return list(pool.map(fn, tasks, chunksize=chunksize))

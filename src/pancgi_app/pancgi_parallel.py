import multiprocessing
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait


def completed_tasks(function, tasks, workers, initializer=None, initargs=()):
    if type(workers) is not int or workers < 1:
        raise ValueError('Worker count must be a positive integer')
    source = iter(enumerate(tasks))
    pool = ProcessPoolExecutor(max_workers=workers,
        mp_context=multiprocessing.get_context('spawn'),
        initializer=initializer, initargs=initargs)
    pending = {}
    try:
        for _ in range(workers):
            item = next(source, None)
            if item is None:
                break
            ordinal, task = item
            pending[pool.submit(function, task)] = ordinal
        while pending:
            done, _ = wait(pending, return_when=FIRST_COMPLETED)
            results = [(pending.pop(future), future.result()) for future in done]
            for ordinal, result in results:
                yield ordinal, result
                item = next(source, None)
                if item is not None:
                    index, task = item
                    pending[pool.submit(function, task)] = index
    finally:
        for future in pending:
            future.cancel()
        pool.shutdown(wait=True, cancel_futures=True)

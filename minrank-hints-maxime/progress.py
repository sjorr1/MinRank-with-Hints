"""Simple, flushed phase messages for terminal experiments."""
from contextlib import contextmanager
from time import perf_counter


@contextmanager
def activity(label, enabled=False):
    """Print a phase start and completion, including its duration in seconds.

    Messages are flushed immediately. While a SageMath operation is running,
    the last message identifies the active phase. An interrupted or failed
    operation prints 'stopped' before propagating the exception to the caller.
    """
    if not enabled:
        yield
        return
    start = perf_counter()
    print(f'{label}...', flush=True)
    outcome = 'finished'
    try:
        yield
    except BaseException:
        outcome = 'stopped'
        raise
    finally:
        print(f'{label}: {outcome} after {perf_counter()-start:.2f} s.', flush=True)

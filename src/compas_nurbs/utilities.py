# TODO: move those functions to compas

def linspace(start, stop, num, endpoint=True):
    """ Returns a list of evenly spaced numbers over a specified interval.

    Inspired from Numpy's linspace function: https://github.com/numpy/numpy/blob/master/numpy/core/function_base.py

    Parameters
    ----------
    start : float
        The starting value of the sequence.
    stop : float
        The end value of the sequence.
    num : int
        Number of samples to generate.
    endpoint : bool, optional
        If True, `stop` is the last sample. Otherwise, it is not included.
        Default is True.

    Returns
    -------
    samples : list
        There are `num` equally spaced samples in the closed interval
        ``[start, stop]`` or the half-open interval ``[start, stop)``
        (depending on whether `endpoint` is True or False).

    Examples
    --------
    >>> linspace(2.0, 3.0, 5)
    [2.0, 2.25, 2.5, 2.75, 3.0]
    >>> linspace(2.0, 3.0, 5, endpoint=False)
    [2.0, 2.2, 2.4, 2.6, 2.8]
    """
    if num < 0:
        raise ValueError("Number of samples, %s, must be non-negative." % num)
    start = float(start)
    stop = float(stop)
    if abs(start - stop) <= 10e-8:
        return [start]
    num = int(num)
    if num > 1:
        div = (num - 1) if endpoint else num
        delta = stop - start
        return [start + (float(x) * float(delta) / float(div)) for x in range(num)]
    else:
        return [start]


def unflatten(array, n):
    if len(array) % n:
        raise ValueError("The length of the array must be a factor of n: %d %% %d == 0" % (len(array), n))
    return [array[i:i + n] for i in range(0, len(array), n)]

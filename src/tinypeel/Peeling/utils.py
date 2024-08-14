from numba import jit
import numpy as np


@jit(nopython=True, nogil=True, fastmath=True)
def summing(vec):
    """Summing the input `vec` over axis 0, equivalent to np.sum(vec, axis=0).

    :param vec: a vector
    :type vec: numpy array with dimension > 1
    :return: a vector
    :rtype: numpy array with one less dimension compared to `vec`
    """
    total = np.zeros(vec.shape[1:], dtype=np.float32)
    for i in range(vec.shape[0]):
        total += vec[i]
    return total


@jit(nopython=True, nogil=True, fastmath=True)
def summing_twice(vec):
    """Summing the input `vec` over axis 0 and 1, equivalent to np.sum(vec, axis=(0, 1)).

    :param vec: a vector
    :type vec: numpy array with dimension > 2
    :return: a vector
    :rtype: numpy array with one less dimension compared to `vec`
    """
    total = np.zeros(vec.shape[2:], dtype=np.float32)
    for i in range(vec.shape[0]):
        for j in range(vec.shape[1]):
            total += vec[i, j]
    return total

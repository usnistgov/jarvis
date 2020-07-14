"""Set of useful utility functions."""

from collections import OrderedDict
from collections import defaultdict
import random


def get_counts(array=['W', 'W', 'Mo', 'Mo', 'S', 'S']):
    """
    Get number of unique elements and their counts.

    Uses OrderedDict.

    Args:
         array of elements

    Returns:
         ordereddict, e.g.OrderedDict([('W', 2), ('Mo', 2), ('S', 2)])
    """
    uniqe_els = []
    for i in array:
        if i not in uniqe_els:
            uniqe_els.append(i)
    info = OrderedDict()
    for i in uniqe_els:
        info.setdefault(i, 0)
    for i in array:
        info[i] += 1
    return info


def gcd(a, b):
    """Calculate the Greatest Common Divisor of a and b.

    Unless b==0, the result will have the same sign as b (so that when
    b is divided by it, the result comes out positive).
    """
    while b:
        a, b = b, a % b
    return a


def ext_gcd(a, b):
    """GCD module from ase."""
    if b == 0:
        return 1, 0
    elif a % b == 0:
        return 0, 1
    else:
        x, y = ext_gcd(b, a % b)
        return y, x - y * (a // b)


def rand_select(x=[]):
    """Select randomly with index info."""
    info = {}
    for i, ii in enumerate(x):
        info.setdefault(ii, []).append(i)
    selected = {}
    for i, j in info.items():
        chosen = random.choice(j)
        selected.setdefault(i, chosen)
    return selected


def rec_dict():
    """Make a recursion dictionary."""
    return defaultdict(rec_dict)

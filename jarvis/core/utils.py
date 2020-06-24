"""
Set of useful utility functions
"""
from collections import defaultdict, OrderedDict
import numpy as np

def get_counts(array=['W','W','Mo','Mo','S','S']):
    """
    Get number of unique elements and their counts
    
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
      info.setdefault(i,0)
    for i in array:
      info[i]+=1
    return info

def gcd(a, b):
    """Calculate the Greatest Common Divisor of a and b.

    Unless b==0, the result will have the same sign as b (so that when
    b is divided by it, the result comes out positive).
    """
    while b:
        a, b = b, a % b
    return a


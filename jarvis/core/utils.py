"""Set of useful utility functions."""

from collections import OrderedDict
from collections import defaultdict
import random
import numpy as np
import math
import xmltodict


def xml_to_dict(fname):
    """Parse XML file."""
    with open(fname, "r") as f:
        data = xmltodict.parse(f.read())
    return data


def get_counts(array=["W", "W", "Mo", "Mo", "S", "S"]):
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


def random_colors(number_of_colors=110):
    """Generate random colors for atom coloring."""
    colors = [
        "#" + "".join([random.choice("0123456789ABCDEF") for j in range(6)])
        for i in range(number_of_colors)
    ]
    color_dict = {}
    for i, ii in enumerate(colors):
        color_dict[i] = ii
    return color_dict


def get_angle(
    a=np.array([1, 2, 3]), b=np.array([4, 5, 6]), c=np.array([7, 8, 9])
):
    """Get angle between three vectors."""
    # theta = argcos(x.y/(|x||y|))
    cos = np.dot((a - b), (c - b)) / (
        np.linalg.norm((a - b)) * np.linalg.norm((c - b))
    )
    if cos <= -1.0:
        cos = cos + 0.000001
    if cos >= 1.0:
        cos = cos - 0.000001
    angle = math.degrees(math.acos(cos))
    return angle


def recast_array_on_uniq_array_elements(
    uniq=["Si", "Al", "O"],
    arr=["Si", "Si", "Al", "Al", "Si", "O", "O", "O", "O"],
):
    """Recast array on uniq array elements."""
    info = {}
    for i, ii in enumerate(uniq):
        for j, jj in enumerate(arr):
            if ii == jj:
                info.setdefault(ii, []).append(j)
    return info


def lorentzian(x, y, x0, gamma):
    """Get Lorentzian of a function."""
    return (y / math.pi) * (
        (0.5 * gamma) / ((x - x0) ** 2 + (0.5 * gamma) ** 2)
    )


# color_dict=random_colors()
def stringdict_to_xml(d={}, enforce_string=False):
    """Convert string dictionary to XML."""
    line = ""
    for i, j in d.items():
        if enforce_string:
            line += "<" + str(i) + ">'" + str(j) + "'</" + str(i) + ">"
        else:
            line += "<" + str(i) + ">" + str(j) + "</" + str(i) + ">"
    return line


def array_to_string(arr=[]):
    """Convert 1D arry to string."""
    return ",".join(map(str, arr))


def chunks(lst, n):
    """Split successive n-sized chunks from list."""
    x = []
    for i in range(0, len(lst), n):
        x.append(lst[i : i + n])
    return x


def check_match(a, b, tol=1e-8):
    """Check if a and b are the same, taking into account PBCs."""
    if abs(a[0] - b[0]) < tol or abs(abs(a[0] - b[0]) - 1) < tol:
        if abs(a[1] - b[1]) < tol or abs(abs(a[1] - b[1]) - 1) < tol:
            if abs(a[2] - b[2]) < tol or abs(abs(a[2] - b[2]) - 1) < tol:
                return True
    return False


# def is_xml_valid(xsd="jarvisdft.xsd", xml="JVASP-1002.xml"):
#    """Check if XML is valid."""
#    xml_file = etree.parse(xml)
#    xml_validator = etree.XMLSchema(file=xsd)
#    is_valid = xml_validator.validate(xml_file)
#    return is_valid

"""Set of useful utility functions."""

from collections import OrderedDict
from collections import defaultdict
import random
import numpy as np
import math
import xmltodict
import re
import requests


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


def update_dict(main={}, extra={}):
    """Return update dictionary."""
    # Helper function for dict.update method
    tmp = main.copy()
    for i, j in extra.items():
        tmp[i] = j
    return tmp


def get_new_coord_for_xyz_sym(frac_coord=[], xyz_string=""):
    """Obtain new coord from xyz string."""
    affine_matrix = parse_xyz_string(xyz_string)
    coord = operate_affine(frac_coord, affine_matrix)
    coord = np.array([i - math.floor(i) for i in coord])
    return coord


def check_duplicate_coords(coords=[], coord=[]):
    """Check if a coordinate exists."""
    positive = False
    for i in coords:
        tmp = check_match(i, coord)
        if tmp:
            positive = True
    return positive


def parse_xyz_string(xyz_string):
    """
    Convert xyz info to translation and rotation vectors.

    Adapted from pymatgen.
    Args:
        xyz_string: string of the form 'x, y, z', '-x, -y, z',
            '-2y+1/2, 3x+1/2, z-y+1/2', etc.
    Returns:
        translation and rotation vectors.
    """
    rot_matrix = np.zeros((3, 3))
    trans = np.zeros(3)
    toks = xyz_string.strip().replace(" ", "").lower().split(",")
    re_rot = re.compile(r"([+-]?)([\d\.]*)/?([\d\.]*)([x-z])")
    re_trans = re.compile(r"([+-]?)([\d\.]+)/?([\d\.]*)(?![x-z])")
    for i, tok in enumerate(toks):
        # build the rotation matrix
        for m in re_rot.finditer(tok):
            factor = -1 if m.group(1) == "-" else 1
            if m.group(2) != "":
                factor *= (
                    float(m.group(2)) / float(m.group(3))
                    if m.group(3) != ""
                    else float(m.group(2))
                )
            j = ord(m.group(4)) - 120
            rot_matrix[i, j] = factor
        # build the translation vector
        for m in re_trans.finditer(tok):
            factor = -1 if m.group(1) == "-" else 1
            num = (
                float(m.group(2)) / float(m.group(3))
                if m.group(3) != ""
                else float(m.group(2))
            )
            trans[i] = num * factor
    affine_matrix = np.eye(4)
    affine_matrix[0:3][:, 0:3] = rot_matrix
    affine_matrix[0:3][:, 3] = trans
    return affine_matrix


def operate_affine(cart_coord=[], affine_matrix=[]):
    """Operate affine method."""
    affine_point = np.array([cart_coord[0], cart_coord[1], cart_coord[2], 1])
    return np.dot(np.array(affine_matrix), affine_point)[0:3]


def gaussian(x, sigma):
    """Get Gaussian profile."""
    return np.exp(-(x ** 2) / (2 * sigma ** 2))


def lorentzian2(x, gamma):
    """Get Lorentziann profile."""
    return (
        gamma
        / 2
        / (np.pi * (x ** 2 + (gamma / 2) ** 2))
        / (2 / (np.pi * gamma))
    )


def digitize_array(values=[], max_len=10):
    """Digitze an array."""
    has_float = False in [float(i).is_integer() for i in values]
    if has_float:
        arr = np.array([float(i) for i in values])
        max_val = max(arr)
        min_val = min(arr)
        bins = np.arange(1, max_len + 1) * (max_val - min_val) / 10
        values = np.digitize(arr, bins)
    return values


def bond_angle(
    dist1,
    dist2,
    bondx1,
    bondx2,
    bondy1,
    bondy2,
    bondz1,
    bondz2,
):
    """Get an angle."""
    nm = dist1 * dist2
    rrx = bondx1 * bondx2
    rry = bondy1 * bondy2
    rrz = bondz1 * bondz2
    cos = (rrx + rry + rrz) / (nm)
    if cos <= -1.0:
        cos = cos + 0.000001
    if cos >= 1.0:
        cos = cos - 0.000001
    deg = math.degrees(math.acos(cos))
    return deg


def check_url_exists(
    url="https://www.ctcms.nist.gov/~knc6/static/JARVIS-DFT/JVASP-77580.xml",
):
    """Check if a url exists."""
    request = requests.get(url)
    if request.status_code == 200:
        return True
    else:
        return False


def volumetric_grid_reshape(data=[], final_grid=[50, 50, 50]):
    """Reshape volumetric data."""
    import torch

    data = torch.tensor(data).unsqueeze(0).unsqueeze(0)
    new_data = (
        torch.nn.functional.interpolate(
            data,
            size=final_grid,
            scale_factor=None,
            mode="trilinear",
            align_corners=True,
            recompute_scale_factor=None,
        )
        .squeeze()
        .squeeze()
    )
    return new_data.numpy()


def cos_formula(a, b, c):
    """Get angle between three edges for oblique triangles."""
    res = (a ** 2 + b ** 2 - c ** 2) / (2 * a * b)
    res = -1.0 if res < -1.0 else res
    res = 1.0 if res > 1.0 else res
    return np.arccos(res)


# def is_xml_valid(xsd="jarvisdft.xsd", xml="JVASP-1002.xml"):
#   """Check if XML is valid."""
#   xml_file = etree.parse(xml)
#   xml_validator = etree.XMLSchema(file=xsd)
#   is_valid = xml_validator.validate(xml_file)
#   return is_valid

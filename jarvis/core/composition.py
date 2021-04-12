"""Modules handling chemical composition."""
# from math import gcd
import string
from jarvis.core.specie import Specie
from collections import OrderedDict
from collections import defaultdict
from jarvis.core.utils import gcd
import re
import numpy as np


class Composition(object):
    """Generate composition python object."""

    def __init__(self, content={}, sort=False):
        """
        Initialize with a dictionary.

        This class is used in Atoms and related objects.

        Args:
            content: dictioary.

            sort: whether to sort alphabetically.

        >>> from composition import Composition
        >>> comp = {"Li": 2, "O": 4}
        >>> cc = Composition(comp)
        >>> weight = round(cc.weight,4)
        >>> print(cc.prototype,cc.formula,cc.reduced_formula,weight)
        AB2 Li2O4 LiO2 77.8796
        """
        if sort:
            content = OrderedDict(
                sorted(content.items(), key=lambda x: (x[0]))
            )
        self._content = content

    @staticmethod
    def from_string(value, sort=True):
        """Generate composition from string."""
        re_formula = re.compile(r"([A-Z][a-z]?)([0-9\.]*)")
        d = defaultdict(float)
        for elt, amt in re_formula.findall(value):
            if elt in ["D", "T"]:
                elt = "H"
            if amt == "":
                d[elt] = 1
            elif float(amt).is_integer():
                d[elt] += int(float(amt))
            else:
                d[elt] += float(amt)

        comp = Composition(dict(d), sort=sort)
        return comp

    def reduce(self):
        """Reduce chemical formula."""
        repeat = 0
        for specie, count in self._content.items():
            if repeat == 0:
                repeat = count
            else:
                repeat = gcd(count, repeat)
        reduced = {}
        for specie, count in self._content.items():
            reduced[specie] = count // repeat
        return reduced, repeat

    @property
    def prototype(self):
        """Get chemical prototypes such as A, AB etc."""
        proto = ""
        all_upper = string.ascii_uppercase

        reduced, repeat = self.reduce()
        N = 0
        for specie, count in reduced.items():
            proto = proto + str(all_upper[N]) + str(round(count, 3))
            N = N + 1
        return proto.replace("1", "")

    def to_dict(self):
        """Return dictionary format."""
        return self._content

    @property
    def nspecies(self):
        """Return number of species."""
        return len(list(self.to_dict().keys()))

    @property
    def search_string(self):
        """Return JARVIS sorted search string."""
        return "-".join(sorted(list(self.to_dict().keys())))

    @classmethod
    def from_dict(self, d={}):
        """Load the class from a dictionary."""
        return Composition(content=d)

    @property
    def reduced_formula(self):
        """Get reduced formula."""
        form = ""
        reduced, repeat = self.reduce()
        X = {}
        for i, j in reduced.items():
            X[i] = Specie(i).X
        Y = dict(sorted(X.items(), key=lambda item: item[1]))
        Z = {}
        for i, j in Y.items():
            Z[i] = reduced[i]
        for specie, count in Z.items():
            if float(count).is_integer():
                form = form + specie + str(int(count))
            else:
                form = form + specie + str(count)

        return form.replace("1", "")

    @property
    def formula(self):
        """Get total chemical formula."""
        form = ""
        for specie, count in self._content.items():
            if float(count).is_integer():
                form = form + str(specie) + str(int(count))
            else:
                form = form + str(specie) + str(count)
        return form.replace("1", "")

    @property
    def atomic_fraction(self):
        """Get atomic fraction."""
        comp_dict = self.to_dict()
        tot = sum(comp_dict.values())
        new_dict = OrderedDict()
        for i, j in comp_dict.items():
            new_dict[i] = j / tot
        return new_dict

    @property
    def atomic_fraction_array(self):
        """Get atomic fraction array."""
        nelements = len(list(Specie()._data.keys()))
        frac_arr = np.zeros(nelements)
        fracs = self.atomic_fraction
        for i, j in fracs.items():
            frac_arr[int(Specie(i).Z) - 1] = j
        return frac_arr

    @property
    def weight(self):
        """Get atomic weight."""
        wt = 0.0
        for specie, count in self._content.items():
            wt = wt + Specie(specie).atomic_mass * count
        return wt

    def __repr__(self):
        """Representation of the class."""
        return str(self._content)


"""
if __name__ == "__main__":
    c=Composition.from_string('Al2O3Al5Co6O1')
    print ('Al2O3Al5Co6O1',c.formula, c.reduced_formula)
    import sys
    sys.exit()
    comp = {"Li": 2, "O": 4}
    cc = Composition(comp)
    print("dict", cc.to_dict())
    x, y = cc.reduce()
    # print(x, y)
    proto = cc.prototype
    print(proto, cc.formula, cc.reduced_formula)
"""

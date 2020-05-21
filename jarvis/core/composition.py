"""
Modules handling chemical composition
"""
# from math import gcd
import string
from jarvis.core.specie import Specie
from collections import OrderedDict

def gcd(a, b):
    """Calculate the Greatest Common Divisor of a and b.

    Unless b==0, the result will have the same sign as b (so that when
    b is divided by it, the result comes out positive).
    """
    while b:
        a, b = b, a % b
    return a


class Composition(object):
    def __init__(self, content={}, sort=False):
        """
        >>> from composition import Composition
        >>> comp = {"Li": 2, "O": 4}
        >>> cc = Composition(comp)
        >>> print(cc.prototype, cc.formula, cc.reduced_formula, round(cc.weight,4))
        AB2 Li2O4 LiO2 77.8796
        """
        if sort:
            content = OrderedDict(sorted(content.items(), key=lambda x: (x[0])))
        self._content = content

    def reduce(self):
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
        proto = ""
        all_upper = string.ascii_uppercase

        reduced, repeat = self.reduce()
        N = 0
        for specie, count in reduced.items():
            proto = proto + str(all_upper[N]) + str(round(count, 3))
            N = N + 1
        return proto.replace("1", "")

    def to_dict(self):
        return self._content

    @property
    def reduced_formula(self):
        form = ""
        reduced, repeat = self.reduce()
        for specie, count in reduced.items():
            form = form + specie + str(count)
        return form.replace("1", "")

    @property
    def formula(self):
        form = ""
        for specie, count in self._content.items():
            form = form + str(specie) + str(count)
        return form.replace("1", "")

    @property
    def weight(self):
        wt = 0.0
        for specie, count in self._content.items():
            wt = wt + Specie(specie).atomic_mass * count
        return wt

    def __repr__(self):
        return str(self._content)


"""
if __name__ == "__main__":
    comp = {"Li": 2, "O": 4}
    cc = Composition(comp)
    print("dict", cc.to_dict())
    x, y = cc.reduce()
    # print(x, y)
    proto = cc.prototype
    print(proto, cc.formula, cc.reduced_formula)
"""

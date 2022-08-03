"""Module to get elemnt fraction based array for a formula."""
# Maximum upto 103 elements
import numpy as np
from jarvis.core.composition import Composition
from jarvis.core.specie import Specie


def get_element_fraction_desc(formula="SiO2", max_nelements=103):
    """Get element fraction."""
    x = np.zeros(max_nelements)
    fracs = Composition.from_string(formula).atomic_fraction
    for i, j in fracs.items():
        # -1 because python array srats from 0
        x[Specie(i).Z - 1] = j
    return x

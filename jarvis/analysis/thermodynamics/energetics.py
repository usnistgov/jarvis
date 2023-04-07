"""Get formation energy, convex hull etc.."""

import os
from scipy.spatial import ConvexHull
import numpy as np
from jarvis.db.figshare import data
from jarvis.core.atoms import Atoms
from jarvis.db.jsonutils import loadjson
from collections import OrderedDict
from jarvis.core.composition import Composition
import re

# import matplotlib.pyplot as plt


def get_optb88vdw_energy():
    """Get OptB88vdW energy per atoms for elements."""
    return loadjson(os.path.join(os.path.dirname(__file__), "unary.json"))


def get_unary_qe_tb_energy():
    """Get elemental chemical potential for GBRV tight-binding project."""
    return loadjson(
        os.path.join(os.path.dirname(__file__), "unary_qe_tb.json")
    )


def isfloat(value):
    """Check if a number is float.

    TODO: replace with isinstance.
    """
    try:
        float(value)
        return True
    except ValueError:
        return False


def unary_energy(el="Na", chem_pots=None):
    """Provide energy per atoms of an element."""
    if chem_pots is None:
        chem_pots = get_optb88vdw_energy()
    en = "na"
    for i, j in chem_pots.items():
        if str(i) == str(el):
            en = j["energy"]
    return en


def form_enp(atoms=None, total_energy=None, chem_pots=None):
    """
    Calculate formation energy given the total energy and the atoms object.

    Currently for OptB88vdW functional based chemical potential implemented
    but can be generalized by changing unary_energy.
    """
    dd = atoms.composition.to_dict()
    # print ('dd',dd)
    ref = 0.0
    for k, v in dd.items():
        e1 = unary_energy(el=k, chem_pots=chem_pots)
        # print (k,v,e1,total_energy)
        if e1 == "na":
            ref = "na"
            ValueError("Element reference does not exist", e1)

        else:
            ref = ref + float(v) * float(e1)
    if isfloat(ref):
        form_en = float(total_energy) - float(ref)
        form_en = round(float(form_en) / float(atoms.num_atoms), 5)
    return form_en


def get_twod_defect_energy(vrun="", jid="", atom=""):
    """Get mono 2D defect formation energy with OptB88vdW data."""
    dft2d = data("dft_2d")

    def get_enp_jid(jid=""):
        for i in dft2d:
            if i["jid"] == jid:
                return (
                    i["optb88vdw_total_energy"]
                    / Atoms.from_dict(i["atoms"]).num_atoms
                )

        # dir='JVASP-667_C_C_c'
        # tmp=dir.split('_')
        # jid=tmp[0]
        # atom=tmp[2]

    strt = vrun.all_structures[-1]
    natoms = strt.num_atoms
    fin_en = vrun.final_energy
    chem_pot = unary_energy(atom)
    bulk_en_pa = get_enp_jid(jid)
    Ef = fin_en - (natoms + 1) * bulk_en_pa + chem_pot
    return Ef


class PhaseDiagram:
    """Module for phase diagram."""

    def __init__(
        self,
        entries,
        verbose=False,
        only_plot_stable=False,
        only_label_stable=False,
    ):
        """Initialize Phase-diagram."""
        # Adapted from ASE
        self.species = OrderedDict()
        # List of formula,formation energy,JID etc.
        self.entries = entries
        self.entries_dict = []
        self.verbose = verbose
        self.only_plot_stable = only_plot_stable
        self.only_label_stable = only_label_stable
        for i in self.entries:
            name = i[0]
            energy = i[1]
            # jid = i[2]
            count = Composition.from_string(name).to_dict()

            natoms = 0
            for symbol, n in count.items():
                natoms += n
                if symbol not in self.species:
                    self.species[symbol] = len(self.species)
            self.entries_dict.append((count, energy, name, natoms))

        ns = len(self.species)
        self.symbols = [None] * ns
        for symbol, id in self.species.items():
            self.symbols[id] = symbol

        if verbose:
            print("Species:", ", ".join(self.symbols))
            print("Entries:", len(self.entries_dict))
            for i, (count, energy, name, natoms) in enumerate(
                self.entries_dict
            ):
                print("{:<5}{:10}{:10.3f}".format(i, name, energy))

        self.points = np.zeros((len(self.entries_dict), ns + 1))
        for s, (count, energy, name, natoms) in enumerate(self.entries_dict):
            for symbol, n in count.items():
                self.points[s, self.species[symbol]] = n / natoms
            self.points[s, -1] = energy  # / natoms

        if len(self.points) == ns:
            # Simple case that qhull would choke on:
            self.simplices = np.arange(ns).reshape((1, ns))
            self.hull = np.ones(ns, bool)
        else:
            # print("self.points[:, 1:]",self.points[:, 1:])
            hull = ConvexHull(self.points[:, 1:])

            # Find relevant simplices:
            ok = hull.equations[:, -2] < 0
            self.simplices = hull.simplices[ok]

            # Create a mask for those points that are on the convex hull:
            self.hull = np.zeros(len(self.points), bool)
            for simplex in self.simplices:
                self.hull[simplex] = True

    def energy_above_hull(self, entry=[]):
        """Find energy above hull."""
        formula = entry[0]
        form_enp = entry[1]
        kwargs = Composition.from_string(formula).to_dict()

        point = np.zeros(len(self.species))
        N = 0
        for symbol, n in kwargs.items():
            point[self.species[symbol]] = n
            N += n
        # print ('N',N)
        # Find coordinates within each simplex:
        X = self.points[self.simplices, 1:-1] - point[1:] / N

        # Find the simplex with positive coordinates that sum to
        # less than one:
        eps = 1e-15
        for i, Y in enumerate(X):
            try:
                x = np.linalg.solve((Y[1:] - Y[:1]).T, -Y[0])
            except np.linalg.linalg.LinAlgError:
                continue
            if (x > -eps).all() and x.sum() < 1 + eps:
                break
        else:
            assert False, X

        indices = self.simplices[i]
        points = self.points[indices]

        scaledcoefs = [1 - x.sum()]
        scaledcoefs.extend(x)
        #         print('scaledcoefs',scaledcoefs)
        #         print('points[:, -1]',points[:, -1])
        energy = np.dot(scaledcoefs, points[:, -1])  # *N

        coefs = []
        results = []
        for coef, s in zip(scaledcoefs, indices):
            count, e, name, natoms = self.entries_dict[s]
            coef *= N / natoms
            coefs.append(coef)
            results.append((name, coef, e))

        #         if self.verbose:
        #             print_results(results)
        e_above_hull = form_enp - energy
        return e_above_hull, energy, indices, np.array(coefs)

    def get_ehull_all(self):
        """Find energy above hull for all entries."""
        info = []
        for i in self.entries:
            # print('ent',i)
            ehull, energy, indices, coefs = self.energy_above_hull(
                entry=[i[0], i[1]]
            )
            info.append([i, ehull])
        return info

    def plot(self, ax=None, dims=None, show=False):
        """Make 2-d or 3-d plot of datapoints and convex hull.

        Default is 2-d for 2- and 3-component diagrams and 3-d for a
        4-component diagram.
        """
        import matplotlib.pyplot as plt

        N = len(self.species)

        if dims is None:
            if N <= 3:
                dims = 2
            else:
                dims = 3

        if ax is None:
            projection = None
            if dims == 3:
                projection = "3d"
                from mpl_toolkits.mplot3d import Axes3D

                Axes3D  # silence pyflakes
            fig = plt.figure()
            ax = fig.add_subplot(projection=projection)
        else:
            if dims == 3 and not hasattr(ax, "set_zlim"):
                raise ValueError(
                    "Cannot make 3d plot unless axes projection " "is 3d"
                )

        if dims == 2:
            if N == 2:
                self.plot2d2(ax)
            elif N == 3:
                self.plot2d3(ax)
            else:
                raise ValueError(
                    "Can only make 2-d plots for 2 and 3 " "component systems!"
                )
        else:
            if N == 3:
                self.plot3d3(ax)
            elif N == 4:
                self.plot3d4(ax)
            else:
                raise ValueError(
                    "Can only make 3-d plots for 3 and 4 " "component systems!"
                )

        if show:
            plt.show()
        return ax

    def plot2d2(self, ax=None):
        """Get 2D plot."""
        x, e = self.points[:, 1:].T
        names = [
            re.sub(r"(\d+)", r"$_{\1}$", ref[2]) for ref in self.entries_dict
        ]
        hull = self.hull
        simplices = self.simplices
        xlabel = self.symbols[1]
        ylabel = "energy [eV/atom]"
        extra = -min(e) / 10
        if ax:
            for i, j in simplices:
                ax.plot(x[[i, j]], e[[i, j]], "-b")
            ax.plot(x[hull], e[hull], "sg")
            if not self.only_plot_stable:
                ax.plot(x[~hull], e[~hull], "or")

            if self.only_plot_stable or self.only_label_stable:
                x = x[self.hull]
                e = e[self.hull]
                names = [name for name, h in zip(names, self.hull) if h]
            for a, b, name in zip(x, e, names):
                if b <= extra:
                    ax.text(a, b, name, ha="center", va="top")

            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)

            ax.set_ylim([min(e) - extra, extra])
        return (x, e, names, hull, simplices, xlabel, ylabel)

    def plot2d3(self, ax=None):
        """Get 2D plot for ternaries."""
        x, y = self.points[:, 1:-1].T.copy()
        x += y / 2
        y *= 3 ** 0.5 / 2
        names = [
            re.sub(r"(\d+)", r"$_{\1}$", ref[2]) for ref in self.entries_dict
        ]
        hull = self.hull
        simplices = self.simplices

        if ax:
            for i, j, k in simplices:
                ax.plot(x[[i, j, k, i]], y[[i, j, k, i]], "-b")
            ax.plot(x[hull], y[hull], "og")
            if not self.only_plot_stable:
                ax.plot(x[~hull], y[~hull], "sr")
            if self.only_plot_stable or self.only_label_stable:
                x = x[self.hull]
                y = y[self.hull]
                names = [name for name, h in zip(names, self.hull) if h]
            for a, b, name in zip(x, y, names):
                ax.text(a, b, name, ha="center", va="top")
        ax.set_xticks([])
        ax.set_yticks([])
        ax.axis("off")
        return (x, y, names, hull, simplices)

    def plot3d3(self, ax):
        """Get 3D plot for ternaries."""
        x, y, e = self.points[:, 1:].T

        ax.scatter(x[self.hull], y[self.hull], e[self.hull], c="g", marker="o")
        if not self.only_plot_stable:
            ax.scatter(
                x[~self.hull], y[~self.hull], e[~self.hull], c="r", marker="s"
            )

        for a, b, c, ref in zip(x, y, e, self.entries_dict):
            name = re.sub(r"(\d+)", r"$_{\1}$", ref[2])
            ax.text(a, b, c, name, ha="center", va="bottom")

        for i, j, k in self.simplices:
            ax.plot(
                x[[i, j, k, i]], y[[i, j, k, i]], zs=e[[i, j, k, i]], c="b"
            )

        ax.set_xlim3d(0, 1)
        ax.set_ylim3d(0, 1)
        ax.view_init(azim=115, elev=30)
        ax.set_xlabel(self.symbols[1])
        ax.set_ylabel(self.symbols[2])
        ax.set_zlabel("energy [eV/atom]")
        ax.set_xticks([])
        ax.set_yticks([])
        ax.axis("off")

    def plot3d4(self, ax):
        """Get 3D plot for quaternaries."""
        x, y, z = self.points[:, 1:-1].T
        a = x / 2 + y + z / 2
        b = 3 ** 0.5 * (x / 2 + y / 6)
        c = (2 / 3) ** 0.5 * z

        ax.scatter(a[self.hull], b[self.hull], c[self.hull], c="g", marker="o")
        if not self.only_plot_stable:
            ax.scatter(
                a[~self.hull], b[~self.hull], c[~self.hull], c="r", marker="s"
            )

        for x, y, z, ref in zip(a, b, c, self.entries_dict):
            name = re.sub(r"(\d+)", r"$_{\1}$", ref[2])
            ax.text(x, y, z, name, ha="center", va="bottom")

        for i, j, k, w in self.simplices:
            ax.plot(
                a[[i, j, k, i, w, k, j, w]],
                b[[i, j, k, i, w, k, j, w]],
                zs=c[[i, j, k, i, w, k, j, w]],
                c="b",
            )

        ax.set_xlim3d(0, 1)
        ax.set_ylim3d(0, 1)
        ax.set_zlim3d(0, 1)
        ax.view_init(azim=115, elev=30)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.axis("off")


def jid_hull(jid="", dataset=[]):
    """Get ehull for a jid and a dataset e.g. dft_3d."""
    from jarvis.db.figshare import data

    if isinstance(dataset, str):
        dataset = data(dataset)
    for i in dataset:
        if i["jid"] == jid:
            system = list(set(i["atoms"]["elements"]))
    z = []
    for i in dataset:
        formula = i["formula"]
        comp = Composition.from_string(formula)
        # atom_frac = comp.atomic_fraction
        all_elms = list(comp.to_dict())
        if (set(all_elms)).issubset(set(system)):
            z.append([i["formula"], i["formation_energy_peratom"], i["jid"]])

    pdj = PhaseDiagram(z)
    # pdj.plot()
    info = pdj.get_ehull_all()
    for i in info:
        if i[0][2] == jid:
            return i


def formula_hull(formula_energy_id=[], dataset=[]):
    """Get ehull for a formula_energy_id pair and a dataset e.g. dft_3d."""
    # e.g. ["Al2O3",-1.0,"JVASP-xyz"]
    # for i in dataset:
    #    if i["jid"] == jid:
    #        system = list(set(i["atoms"]["elements"]))
    from jarvis.db.figshare import data

    if isinstance(dataset, str):
        dataset = data(dataset)

    c = Composition.from_string(formula_energy_id[0])
    system = list(c.to_dict().keys())
    z = []
    z.append(formula_energy_id)

    for i in dataset:
        formula = i["formula"]
        comp = Composition.from_string(formula)
        # atom_frac = comp.atomic_fraction
        all_elms = list(comp.to_dict())
        if (set(all_elms)).issubset(set(system)):
            z.append([i["formula"], i["formation_energy_peratom"], i["jid"]])

    pdj = PhaseDiagram(z)
    # pdj.plot()
    info = pdj.get_ehull_all()
    for i in info:
        if i[0][2] == formula_energy_id[-1]:
            return i

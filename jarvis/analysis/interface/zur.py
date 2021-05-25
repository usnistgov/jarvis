"""Design interface using Zur algorithm and Anderson rule."""

from itertools import product
import numpy as np
from jarvis.core.atoms import add_atoms, fix_pbc
from jarvis.core.lattice import Lattice


class ZSLGenerator(object):
    """
    Uses Zur algorithm to find best matched interfaces.

    This class is modified from pymatgen.
    """

    def __init__(
        self,
        max_area_ratio_tol=0.09,
        max_area=400,
        max_length_tol=0.03,
        max_angle_tol=0.01,
    ):
        """
        Intialize for a specific film and substrate.

        Parameters for the class.
        Args:
            max_area_ratio_tol(float): Max tolerance on ratio of
                super-lattices to consider equal

            max_area(float): max super lattice area to generate in search

            max_length_tol: maximum length tolerance in checking if two
                vectors are of nearly the same length

            max_angle_tol: maximum angle tolerance in checking of two sets
                of vectors have nearly the same angle between them
        """
        self.max_area_ratio_tol = max_area_ratio_tol
        self.max_area = max_area
        self.max_length_tol = max_length_tol
        self.max_angle_tol = max_angle_tol

    def is_same_vectors(self, vec_set1, vec_set2):
        """
        Check two sets of vectors are the same.

        Args:
            vec_set1(array[array]): an array of two vectors

            vec_set2(array[array]): second array of two vectors
        """
        if (
            np.absolute(rel_strain(vec_set1[0], vec_set2[0]))
            > self.max_length_tol
        ):
            return False
        elif (
            np.absolute(rel_strain(vec_set1[1], vec_set2[1]))
            > self.max_length_tol
        ):
            return False
        elif np.absolute(rel_angle(vec_set1, vec_set2)) > self.max_angle_tol:
            return False
        else:
            return True

    def generate_sl_transformation_sets(self, film_area, substrate_area):
        """Generate transformation sets for film/substrate.

        The transformation sets map the film and substrate unit cells to super
        lattices with a maximum area.

        Args:

            film_area(int): the unit cell area for the film.

            substrate_area(int): the unit cell area for the substrate.

        Returns:
            transformation_sets: a set of transformation_sets defined as:
                1.) the transformation matricies for the film to create a
                super lattice of area i*film area
                2.) the tranformation matricies for the substrate to create
                a super lattice of area j*film area
        """
        transformation_indicies = [
            (i, j)
            for i in range(1, int(self.max_area / film_area))
            for j in range(1, int(self.max_area / substrate_area))
            if np.absolute(film_area / substrate_area - float(j) / i)
            < self.max_area_ratio_tol
        ]

        # Sort sets by the square of the matching area and yield in order
        # from smallest to largest
        for x in sorted(transformation_indicies, key=lambda x: x[0] * x[1]):
            yield (
                gen_sl_transform_matricies(x[0]),
                gen_sl_transform_matricies(x[1]),
            )

    def get_equiv_transformations(
        self, transformation_sets, film_vectors, substrate_vectors
    ):
        """
        Apply the transformation_sets to the film and substrate vectors.

        Generate super-lattices and checks if they matches.
        Returns all matching vectors sets.

        Args:
            transformation_sets(array): an array of transformation sets:
                each transformation set is an array with the (i,j)
                indicating the area multipes of the film and subtrate it
                corresponds to, an array with all possible transformations
                for the film area multiple i and another array for the
                substrate area multiple j.

            film_vectors(array): film vectors to generate super lattices.

            substrate_vectors(array): substrate vectors to generate super
                lattices
        """
        for (
            film_transformations,
            substrate_transformations,
        ) in transformation_sets:
            # Apply transformations and reduce using Zur reduce methodology
            films = [
                reduce_vectors(*np.dot(f, film_vectors))
                for f in film_transformations
            ]

            substrates = [
                reduce_vectors(*np.dot(s, substrate_vectors))
                for s in substrate_transformations
            ]

            # Check if equivalant super lattices
            for (f_trans, s_trans), (f, s) in zip(
                product(film_transformations, substrate_transformations),
                product(films, substrates),
            ):
                if self.is_same_vectors(f, s):
                    yield [f, s, f_trans, s_trans]

    def __call__(self, film_vectors, substrate_vectors, lowest=False):
        """Run the ZSL algorithm to generate all possible matching."""
        film_area = vec_area(*film_vectors)
        substrate_area = vec_area(*substrate_vectors)

        # Generate all super lattice comnbinations for a given set of miller
        # indicies
        transformation_sets = self.generate_sl_transformation_sets(
            film_area, substrate_area
        )

        # Check each super-lattice pair to see if they match
        for match in self.get_equiv_transformations(
            transformation_sets, film_vectors, substrate_vectors
        ):
            # Yield the match area, the miller indicies,
            yield self.match_as_dict(
                match[0],
                match[1],
                film_vectors,
                substrate_vectors,
                vec_area(*match[0]),
                match[2],
                match[3],
            )

            if lowest:
                break

    def match_as_dict(
        self,
        film_sl_vectors,
        substrate_sl_vectors,
        film_vectors,
        substrate_vectors,
        match_area,
        film_transformation,
        substrate_transformation,
    ):
        """
        Return dict which contains ZSL match.

        Args:
            film_miller(array)

            substrate_miller(array)
        """
        d = {}
        d["film_sl_vecs"] = np.asarray(film_sl_vectors)
        d["sub_sl_vecs"] = np.asarray(substrate_sl_vectors)
        d["match_area"] = match_area
        d["film_vecs"] = np.asarray(film_vectors)
        d["sub_vecs"] = np.asarray(substrate_vectors)
        d["film_transformation"] = np.asarray(film_transformation)
        d["substrate_transformation"] = np.asarray(substrate_transformation)

        return d


def gen_sl_transform_matricies(area_multiple):
    """
    Generate the transformation matricies.

    Convert a set of 2D vectors into a super
    lattice of integer area multiple as proven
    in Cassels:
    Cassels, John William Scott. An introduction to the geometry of
    numbers. Springer Science & Business Media, 2012.

    Args:
        area_multiple(int): integer multiple of unit cell area for super
        lattice area.

    Returns:
        matrix_list: transformation matricies to covert unit vectors to
        super lattice vectors.
    """
    return [
        np.array(((i, j), (0, area_multiple / i)))
        for i in get_factors(area_multiple)
        for j in range(area_multiple // i)
    ]


def rel_strain(vec1, vec2):
    """Calculate relative strain between two vectors."""
    return fast_norm(vec2) / fast_norm(vec1) - 1


def rel_angle(vec_set1, vec_set2):
    """
    Calculate the relative angle between two vector sets.

    Args:
        vec_set1(array[array]): an array of two vectors.

        vec_set2(array[array]): second array of two vectors.
    """
    return (
        vec_angle(vec_set2[0], vec_set2[1])
        / vec_angle(vec_set1[0], vec_set1[1])
        - 1
    )


def fast_norm(a):
    """Much faster variant of numpy linalg norm."""
    return np.sqrt(np.dot(a, a))


def vec_angle(a, b):
    """Calculate angle between two vectors."""
    cosang = np.dot(a, b)
    sinang = fast_norm(np.cross(a, b))
    return np.arctan2(sinang, cosang)


def vec_area(a, b):
    """Area of lattice plane defined by two vectors."""
    return fast_norm(np.cross(a, b))


def reduce_vectors(a, b):
    """Generate independent and unique basis vectors based on Zur et al."""
    if np.dot(a, b) < 0:
        return reduce_vectors(a, -b)
    if fast_norm(a) > fast_norm(b):
        return reduce_vectors(b, a)
    if fast_norm(b) > fast_norm(np.add(b, a)):
        return reduce_vectors(a, np.add(b, a))
    if fast_norm(b) > fast_norm(np.subtract(b, a)):
        return reduce_vectors(a, np.subtract(b, a))
    return [a, b]


def get_factors(n):
    """Generate all factors of n."""
    for x in range(1, n + 1):
        if n % x == 0:
            yield x


def make_interface(
    film="",
    subs="",
    atol=1,
    ltol=0.05,
    max_area=500,
    max_area_ratio_tol=1.00,
    seperation=3.0,
    vacuum=8.0,
    apply_strain=False,
):
    """
    Use as main function for making interfaces/heterostructures.

    Return mismatch and other information as info dict.

    Args:
       film: top/film material.

       subs: substrate/bottom/fixed material.

       seperation: minimum seperation between two.

       vacuum: vacuum will be added on both sides.
       So 2*vacuum will be added.
    """
    z = ZSLGenerator(
        max_area_ratio_tol=max_area_ratio_tol,
        max_area=max_area,
        max_length_tol=ltol,
        max_angle_tol=atol,
    )
    film = fix_pbc(film.center_around_origin([0, 0, 0]))
    subs = fix_pbc(subs.center_around_origin([0, 0, 0]))
    matches = list(z(film.lattice_mat[:2], subs.lattice_mat[:2], lowest=True))
    info = {}
    info["mismatch_u"] = "na"
    info["mismatch_v"] = "na"
    info["mismatch_angle"] = "na"
    info["area1"] = "na"
    info["area2"] = "na"
    info["film_sl"] = "na"
    info["matches"] = matches
    info["subs_sl"] = "na"
    uv1 = matches[0]["sub_sl_vecs"]
    uv2 = matches[0]["film_sl_vecs"]
    u = np.array(uv1)
    v = np.array(uv2)
    a1 = u[0]
    a2 = u[1]
    b1 = v[0]
    b2 = v[1]
    mismatch_u = np.linalg.norm(b1) / np.linalg.norm(a1) - 1
    mismatch_v = np.linalg.norm(b2) / np.linalg.norm(a2) - 1
    angle1 = (
        np.arccos(np.dot(a1, a2) / np.linalg.norm(a1) / np.linalg.norm(a2))
        * 180
        / np.pi
    )
    angle2 = (
        np.arccos(np.dot(b1, b2) / np.linalg.norm(b1) / np.linalg.norm(b2))
        * 180
        / np.pi
    )
    mismatch_angle = abs(angle1 - angle2)
    area1 = np.linalg.norm(np.cross(a1, a2))
    area2 = np.linalg.norm(np.cross(b1, b2))
    uv_substrate = uv1
    uv_film = uv2
    substrate_latt = Lattice(
        np.array(
            [uv_substrate[0][:], uv_substrate[1][:], subs.lattice_mat[2, :]]
        )
    )
    _, __, scell = subs.lattice.find_matches(
        substrate_latt, ltol=ltol, atol=atol
    )
    film_latt = Lattice(
        np.array([uv_film[0][:], uv_film[1][:], film.lattice_mat[2, :]])
    )
    scell[2] = np.array([0, 0, 1])
    scell_subs = scell
    _, __, scell = film.lattice.find_matches(film_latt, ltol=ltol, atol=atol)
    scell[2] = np.array([0, 0, 1])
    scell_film = scell
    film_scell = film.make_supercell_matrix(scell_film)
    subs_scell = subs.make_supercell_matrix(scell_subs)
    info["mismatch_u"] = mismatch_u
    info["mismatch_v"] = mismatch_v
    print("mismatch_u,mismatch_v", mismatch_u, mismatch_v)
    info["mismatch_angle"] = mismatch_angle
    info["area1"] = area1
    info["area2"] = area2
    info["film_sl"] = film_scell
    info["subs_sl"] = subs_scell
    substrate_top_z = max(np.array(subs_scell.cart_coords)[:, 2])
    substrate_bot_z = min(np.array(subs_scell.cart_coords)[:, 2])
    film_top_z = max(np.array(film_scell.cart_coords)[:, 2])
    film_bottom_z = min(np.array(film_scell.cart_coords)[:, 2])
    thickness_sub = abs(substrate_top_z - substrate_bot_z)
    thickness_film = abs(film_top_z - film_bottom_z)
    sub_z = (
        (vacuum + substrate_top_z)
        * np.array(subs_scell.lattice_mat[2, :])
        / np.linalg.norm(subs_scell.lattice_mat[2, :])
    )
    shift_normal = (
        sub_z / np.linalg.norm(sub_z) * seperation / np.linalg.norm(sub_z)
    )
    tmp = (
        thickness_film / 2 + seperation + thickness_sub / 2
    ) / np.linalg.norm(subs_scell.lattice_mat[2, :])
    shift_normal = (
        tmp
        * np.array(subs_scell.lattice_mat[2, :])
        / np.linalg.norm(subs_scell.lattice_mat[2, :])
    )
    interface = add_atoms(
        film_scell, subs_scell, shift_normal, apply_strain=apply_strain
    ).center_around_origin([0, 0, 0.5])
    combined = interface.center(vacuum=vacuum).center_around_origin(
        [0, 0, 0.5]
    )
    info["interface"] = combined
    return info


def get_hetero_type(A={}, B={}):
    """Provide heterojunction classification using Anderson rule."""
    stack = "na"
    int_type = "na"
    try:
        if A["scf_vbm"] - A["avg_max"] < B["scf_vbm"] - B["avg_max"]:
            stack = "BA"
        else:
            C = A
            D = B
            A = D
            B = C
            stack = "AB"
        vbm_a = A["scf_vbm"] - A["avg_max"]
        vbm_b = B["scf_vbm"] - B["avg_max"]
        cbm_a = A["scf_cbm"] - A["avg_max"]
        cbm_b = B["scf_cbm"] - B["avg_max"]
        if vbm_a < vbm_b and vbm_b < cbm_b and cbm_b < cbm_a:
            int_type = "I"
        elif vbm_a < vbm_b and vbm_b < cbm_a and cbm_a < cbm_b:
            int_type = "II"
        elif vbm_a < cbm_a and cbm_a < vbm_b and vbm_b < cbm_b:
            int_type = "III"
    except Exception:
        pass
    return int_type, stack


def mismatch_strts():
    """
    Return mismatch and other information as info dict.

    Deprecated, preserved here so legacy imports work.
    """
    pass


def get_hetero():
    """Generate heterostructure. Deprecated also."""
    pass


"""
if __name__ == "__main__":
    s1 = Poscar.from_file(
        "MAIN-RELAX-Surf-mp-1821/POSCAR"
    )
    s2 = Poscar.from_file(
        "MAIN-RELAX-Surf-mp-2815/POSCAR"
    )
    print(s1)
    print(s2)
    info = mismatch_strts(film=s1.atoms, subs=s2.atoms)
    # print (info)
    print((s1.__dict__["atoms"].lattice_mat))
    d = s1.__dict__["atoms"]
    f = open("test.json", "wb")
    joblib.dump(s1, f)
    f.close()
    print(get_hetero(s1.atoms, s2.atoms))
    ff = open("test.json", "rb")
    dd = joblib.load(ff)
    ff.close()
    # print (dd.atoms.lattice_mat)
"""

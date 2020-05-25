from itertools import product
import numpy as np
from jarvis.core.atoms import Atoms
from jarvis.io.vasp.inputs import Poscar
from jarvis.core.lattice import Lattice
import json
import joblib, pickle


class ZSLGenerator(object):
    """
    Uses Zur algorithm to find best matched interfaces
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
        Intialize a Zur Super Lattice Generator for a specific film and
            substrate
            
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
        Determine if two sets of vectors are the same within length and angle
        tolerances
        
        Args:
            vec_set1(array[array]): an array of two vectors
            
            vec_set2(array[array]): second array of two vectors
            
        """
        if np.absolute(rel_strain(vec_set1[0], vec_set2[0])) > self.max_length_tol:
            return False
        elif np.absolute(rel_strain(vec_set1[1], vec_set2[1])) > self.max_length_tol:
            return False
        elif np.absolute(rel_angle(vec_set1, vec_set2)) > self.max_angle_tol:
            return False
        else:
            return True

    def generate_sl_transformation_sets(self, film_area, substrate_area):
        """
        Generates transformation sets for film/substrate pair given the
        area of the unit cell area for the film and substrate. The
        transformation sets map the film and substrate unit cells to super
        lattices with a maximum area
        
        Args:
        
            film_area(int): the unit cell area for the film
            
            substrate_area(int): the unit cell area for the substrate
            
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
            yield (gen_sl_transform_matricies(x[0]), gen_sl_transform_matricies(x[1]))

    def get_equiv_transformations(
        self, transformation_sets, film_vectors, substrate_vectors
    ):
        """
        Applies the transformation_sets to the film and substrate vectors
        to generate super-lattices and checks if they matches.
        Returns all matching vectors sets.
        
        Args:
        
            transformation_sets(array): an array of transformation sets:
                each transformation set is an array with the (i,j)
                indicating the area multipes of the film and subtrate it
                corresponds to, an array with all possible transformations
                for the film area multiple i and another array for the
                substrate area multiple j.
                
            film_vectors(array): film vectors to generate super lattices
            
            substrate_vectors(array): substrate vectors to generate super
                lattices
        """

        for (film_transformations, substrate_transformations) in transformation_sets:
            # Apply transformations and reduce using Zur reduce methodology
            films = [
                reduce_vectors(*np.dot(f, film_vectors)) for f in film_transformations
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
        """
        Runs the ZSL algorithm to generate all possible matching
        :return:
        """

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

            # Just want lowest match per direction
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
        Returns dict which contains ZSL match
        
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
    Generates the transformation matricies that convert a set of 2D
    vectors into a super lattice of integer area multiple as proven
    in Cassels:
    Cassels, John William Scott. An introduction to the geometry of
    numbers. Springer Science & Business Media, 2012.
    
    Args:
    
        area_multiple(int): integer multiple of unit cell area for super
        lattice area
        
    Returns:
        matrix_list: transformation matricies to covert unit vectors to
        super lattice vectors
        
    """
    return [
        np.array(((i, j), (0, area_multiple / i)))
        for i in get_factors(area_multiple)
        for j in range(area_multiple // i)
    ]


def rel_strain(vec1, vec2):
    """
    Calculate relative strain between two vectors
    """
    return fast_norm(vec2) / fast_norm(vec1) - 1


def rel_angle(vec_set1, vec_set2):
    """
    Calculate the relative angle between two vector sets
    
    Args:
        vec_set1(array[array]): an array of two vectors
        
        vec_set2(array[array]): second array of two vectors
    """
    return vec_angle(vec_set2[0], vec_set2[1]) / vec_angle(vec_set1[0], vec_set1[1]) - 1


def fast_norm(a):
    """
    Much faster variant of numpy linalg norm
    """
    return np.sqrt(np.dot(a, a))


def vec_angle(a, b):
    """
    Calculate angle between two vectors
    """
    cosang = np.dot(a, b)
    sinang = fast_norm(np.cross(a, b))
    return np.arctan2(sinang, cosang)


def vec_area(a, b):
    """
    Area of lattice plane defined by two vectors
    """
    return fast_norm(np.cross(a, b))


def reduce_vectors(a, b):
    """
    Generate independent and unique basis vectors based on the
    methodology of Zur and McGill
    """
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
    """
    Generate all factors of n
    """
    for x in range(1, n + 1):
        if n % x == 0:
            yield x

  
def make_interface(film='',subs='',atol=1,ltol=0.05):
    """
    Returns mismatch and other information as info dict
    """


    z = ZSLGenerator()

    matches = list(z(film.lattice_mat[:2], subs.lattice_mat[:2], lowest=True))

    info = {}
    info["mismatch_u"] = "na"
    info["mismatch_v"] = "na"
    info["mismatch_angle"] = "na"
    info["area1"] = "na"
    info["area2"] = "na"
    info["film_sl"] = "na"#film
    info['matches']= matches
    info["subs_sl"] = "na"#subs

    uv1 = matches[0]["sub_sl_vecs"]
    uv2 = matches[0]["film_sl_vecs"]
    # uv1=[[-8.52917200e+00,  9.12745800e+00,  3.66344517e-17],[1.27937580e+01, 9.12745800e+00, 1.34228735e-15]]
    # uv2=[[7.02403800e+00, 1.05360570e+01, 1.07524571e-15],[-1.40480760e+01,  7.02403800e+00, -4.30098283e-16]]
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
        np.array([uv_substrate[0][:], uv_substrate[1][:], subs.lattice_mat[2, :]])
    )

    _, __, scell = subs.lattice.find_matches(substrate_latt, ltol=ltol, atol=atol)


    film_latt = Lattice(
        np.array([uv_film[0][:], uv_film[1][:], film.lattice_mat[2, :]])
    )


    scell[2] = np.array([0, 0, 1])
    scell_subs = scell

    _, __, scell = film.lattice.find_matches(film_latt, ltol=ltol, atol=atol)
    scell[2] = np.array([0, 0, 1])
    scell_film = scell


    film = film.make_supercell(scell_film)
    subs = subs.make_supercell(scell_subs)

    lmap = Lattice(
        np.array(
            [subs.lattice_mat[0, :], subs.lattice_mat[1, :], film.lattice_mat[2, :]]
        )
    )

    film = film.center_around_origin()
    subs = subs.center_around_origin()
    info["mismatch_u"] = mismatch_u
    info["mismatch_v"] = mismatch_v
    info["mismatch_angle"] = mismatch_angle
    info["area1"] = area1
    info["area2"] = area2
    info["film_sl"] = film
    info["subs_sl"] = subs

    
    seperation=3
    coords_uniq_sub = np.array(subs.cart_coords)


    coords_uniq_film = np.array(film.cart_coords)

    substrate_top_z = max(np.array(subs.frac_coords)[:, 2])
    substrate_bot_z = min(np.array(subs.frac_coords)[:, 2])

    film_bottom = min(np.array(film.frac_coords)[:, 2])
    film_top = max(np.array(film.frac_coords)[:, 2])

    sub_z = subs.lattice_mat[2, :]
    origin = np.array([0, 0, substrate_top_z])
    shift_normal = sub_z / np.linalg.norm(sub_z) * seperation/ np.linalg.norm(sub_z)

    thickness_sub = abs(substrate_top_z-substrate_bot_z)
    thickness_film = abs(film_top-film_bottom)


    new_coords = []
    lattice_mat = subs.lattice_mat
    elements = []
    for i in subs.frac_coords:
        new_coords.append(i)
    for i in subs.elements:
        elements.append(i)

    for i in film.elements:
        elements.append(i)
    for i in film.frac_coords:
        tmp = i
        tmp[2] = i[2] +shift_normal[2]+(thickness_sub+thickness_film)/2.0


        new_coords.append(tmp)


    interface = Atoms(
        lattice_mat=lattice_mat, elements=elements, coords=new_coords, cartesian=False
    )
    info['interface']=interface
    return info




def mismatch_strts(film=[], subs=[],ltol=0.05, atol=10):
    """
    Returns mismatch and other information as info dict,
    Deprecated, this module will be removes in the next version
    """
    z = ZSLGenerator()
    matches = list(z(film.lattice_mat[:2], subs.lattice_mat[:2], lowest=True))
    info = {}
    info["mismatch_u"] = "na"
    info["mismatch_v"] = "na"
    info["mismatch_angle"] = "na"
    info["area1"] = "na"
    info["area2"] = "na"
    info["film_sl"] = film
    info["subs_sl"] = subs

    uv1 = matches[0]["sub_sl_vecs"]
    uv2 = matches[0]["film_sl_vecs"]
    # uv1=[[-8.52917200e+00,  9.12745800e+00,  3.66344517e-17],[1.27937580e+01, 9.12745800e+00, 1.34228735e-15]]
    # uv2=[[7.02403800e+00, 1.05360570e+01, 1.07524571e-15],[-1.40480760e+01,  7.02403800e+00, -4.30098283e-16]]
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
    uv_mat2d = uv2
    substrate_latt = Lattice(
        np.array([uv_substrate[0][:], uv_substrate[1][:], subs.lattice_mat[2, :]])
    )
    # to avoid numerical issues with find_mapping
    mat2d_fake_c = film.lattice_mat[2, :] / np.linalg.norm(film.lattice_mat[2, :]) * 5.0
    mat2d_latt = Lattice(np.array([uv_mat2d[0][:], uv_mat2d[1][:], mat2d_fake_c]))
    mat2d_latt_fake = Lattice(
        np.array([film.lattice_mat[0, :], film.lattice_mat[1, :], mat2d_fake_c])
    )
    _, __, scell = subs.lattice.find_matches(substrate_latt, ltol=ltol, atol=atol)
    scell[2] = np.array([0, 0, 1])
    # print("subs lattice", subs.lattice_mat)
    # print("scell", scell)
    subs = subs.make_supercell_matrix(scell)
    _, __, scell = mat2d_latt_fake.find_matches(mat2d_latt, ltol=0.05, atol=1)
    scell[2] = np.array([0, 0, 1])
    film = film.make_supercell_matrix(scell)
    # modify the substrate lattice so that the 2d material can be
    # grafted on top of it
    #lmap = Lattice(
    #    np.array(
    #        [subs.lattice_mat[0, :], subs.lattice_mat[1, :], film.lattice_mat[2, :]]
    #    )
    #)
    #film.lattice = lmap

    lmap = Lattice(
        np.array(
            [subs.lattice_mat[0, :], subs.lattice_mat[1, :], film.lattice_mat[2, :]]
        )
    )
    film.lattice = lmap



    # film.modify_lattice(lmap)
    # print ("film",film)
    # print ("subs",subs)

    info["mismatch_u"] = mismatch_u
    info["mismatch_v"] = mismatch_v
    info["mismatch_angle"] = mismatch_angle
    info["area1"] = area1
    info["area2"] = area2
    info["film_sl"] = film
    info["subs_sl"] = subs
    return info


def get_hetero_type(A={}, B={}):
    stack = "na"
    int_type = "na"
    try:
        # if A['phi']>B['phi']:
        if A["scf_vbm"] - A["avg_max"] < B["scf_vbm"] - B["avg_max"]:
            stack = "BA"
        else:
            C = A
            D = B
            A = D
            B = C
            stack = "AB"
            # tmp=B
            # B=A
            # A=tmp
        vbm_a = A["scf_vbm"] - A["avg_max"]
        vbm_b = B["scf_vbm"] - B["avg_max"]
        cbm_a = A["scf_cbm"] - A["avg_max"]
        cbm_b = B["scf_cbm"] - B["avg_max"]
        #  print ('vbm_a,vbm_b,cbm_b,cbm_a',vbm_a,vbm_b,cbm_b,cbm_a)
        if vbm_a < vbm_b and vbm_b < cbm_b and cbm_b < cbm_a:
            int_type = "I"
        elif vbm_a < vbm_b and vbm_b < cbm_a and cbm_a < cbm_b:
            int_type = "II"
        elif vbm_a < cbm_a and cbm_a < vbm_b and vbm_b < cbm_b:
            int_type = "III"
    except:
        pass
    return int_type, stack


def get_hetero(film, substrate, seperation=3.0):
    # unique site coordinates in the substrate top layers
    coords_uniq_sub = np.array(substrate.cart_coords)

    # unique site coordinates in the 2D material bottom layers
    coords_uniq_film = np.array(film.cart_coords)
   
    substrate_top_z = max(np.array(substrate.cart_coords)[:, 2])
    substrate_bot_z = min(np.array(substrate.cart_coords)[:, 2])
    #print ('substrate_top_z',substrate_top_z- substrate_bot_z)
    film_bottom = min(np.array(film.cart_coords)[:, 2])
    film_top = max(np.array(film.cart_coords)[:, 2])
    #print ('film_bottom',film_top-film_bottom)
    # shift normal to the surface by 'seperation'
    sub_z = substrate.lattice_mat[2, :]
    origin = np.array([0, 0, substrate_top_z])
    shift_normal = sub_z / np.linalg.norm(sub_z) * seperation
    #print ('origin', origin)
    #print ('shift_normal', shift_normal)
    thickness_sub = abs(substrate_top_z-substrate_bot_z)
    thickness_film = abs(film_top-film_bottom)
    #shift_normal =  seperation
    #shift_normal = sub_z / np.linalg.norm(sub_z) * seperation
    # generate all possible interfaces, one for each combination of
    # unique substrate and unique 2d materials site in the layers .i.e
    # an interface structure for each parallel shift
    # interface = 2D material + substrate

    new_coords = []
    lattice_mat = substrate.lattice_mat
    elements = []
    for i in substrate.cart_coords:
        new_coords.append(i)
    for i in substrate.elements:
        elements.append(i)

    for i in film.elements:
        elements.append(i)
    for i in film.cart_coords:
        tmp = i
        #tmp[2] = i[2] + (thickness_sub+thickness_film)
        #tmp[2] = i[2] +(thickness_sub+thickness_film)+seperation
        #tmp = tmp  + shift_normal
        tmp = tmp + origin + shift_normal
        #new_coords.append(i)
        new_coords.append(tmp)

    interface = Atoms(
        lattice_mat=lattice_mat, elements=elements, coords=new_coords, cartesian=True
    )

    return interface


"""
if __name__ == "__main__":
    s1 = Poscar.from_file(
        "/rk2/knc6/JARVIS-DFT/2D-1L/POSCAR-mp-1821-1L.vasp_PBEBO/MAIN-RELAX-Surf-mp-1821/POSCAR"
    )
    s2 = Poscar.from_file(
        "/rk2/knc6/JARVIS-DFT/2D-1L/POSCAR-mp-2815-1L.vasp_PBEBO/MAIN-RELAX-Surf-mp-2815/POSCAR"
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

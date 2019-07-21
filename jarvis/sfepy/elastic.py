"""Solve the linear elaticity problem to generate data.

Use SfePy to solve a linear strain problem in 2D with a varying
microstructure on a rectangular grid. The rectangle (cube) is held at
the negative edge (plane) and displaced by 1 on the positive x edge
(plane). Periodic boundary conditions are applied to the other
boundaries.

The microstructure is of shape (n_samples, n_x, n_y) or (n_samples,
n_x, n_y, n_z).


"""

import pytest
import numpy as np
from toolz.curried import pipe, do, first, merge, curry, compose
from toolz.curried import map as map_
from toolz.sandbox.parallel import fold

try:
    import sfepy  # pylint: disable=unused-import; # noqa: F401
except ImportError:  # pragma: no cover
    pytest.importorskip("sfepy")
    raise

from sfepy.base.goptions import goptions
from sfepy.discrete.fem import Field

try:
    from sfepy.discrete.fem import FEDomain as Domain
except ImportError:  # pragma: no cover
    from sfepy.discrete.fem import Domain
from sfepy.discrete import (
    FieldVariable,
    Material,
    Integral,
    Function,
    Equation,
    Equations,
    Problem,
)
from sfepy.terms import Term
from sfepy.discrete.conditions import Conditions, EssentialBC, PeriodicBC
from sfepy.solvers.ls import ScipyDirect
from sfepy.solvers.nls import Newton
import sfepy.discrete.fem.periodic as per
from sfepy.discrete import Functions
from sfepy.mesh.mesh_generators import gen_block_mesh
from sfepy.mechanics.matcoefs import ElasticConstants
from sfepy.base.base import output
from sfepy.discrete.conditions import LinearCombinationBC
from sfepy.mechanics.matcoefs import stiffness_from_lame


goptions["verbose"] = False
output.set_output(quiet=True)


def sequence(*args):
    """Compose functions in order

    Args:
      args: the functions to compose

    Returns:
      composed functions

    """
    return compose(*args[::-1])


@curry
def solve(x_data, elastic_modulus, poissons_ratio, macro_strain=1.0, delta_x=1.0):
    """Solve the elasticity problem

    Args:
      x_data: microstructure with shape (n_samples, n_x, ...)
      elastic_modulus: the elastic modulus in each phase
      poissons_ration: the poissons ratio for each phase
      macro_strain: the macro strain
      delta_x: the grid spacing

    Returns:
      a dictionary of strain, displacement and stress with stress and
      strain of shape (n_samples, n_x, ..., 3) and displacement shape
      of (n_samples, n_x + 1, ..., 2)

    """

    def solve_one_sample(property_array):
        return pipe(
            get_fields(property_array.shape[:-1], delta_x),
            lambda x: get_problem(x, property_array, delta_x, macro_strain),
            lambda x: (x, x.solve()),
            lambda x: get_data(property_array.shape[:-1], *x),
        )

    convert = lambda x: _convert_properties(
        len(x.shape) - 1, elastic_modulus, poissons_ratio
    )[x]

    solve_multiple_samples = sequence(
        do(_check(len(elastic_modulus), len(poissons_ratio))),
        convert,
        map_(solve_one_sample),
        lambda x: zip(*x),
        map_(np.array),
        lambda x: zip(("strain", "displacement", "stress"), tuple(x)),
        dict,
    )

    return solve_multiple_samples(x_data)


def _convert_properties(dim, elastic_modulus, poissons_ratio):
    """Convert from elastic modulus and Poisson's ratio to the Lame
    parameter and shear modulus

    Args:
      dim: whether 2D or 3D
      elastic_modulus: the elastic modulus in each phase
      poissons_ration: the poissons ratio for each phase

    Returns:
      array of shape (n_phases, 2) where for example [1, 0], gives the
      Lame parameter in the second phase


    """
    return pipe(
        zip(elastic_modulus, poissons_ratio),
        map_(
            lambda x: pipe(
                ElasticConstants(young=x[0], poisson=x[1]),
                lambda y: (y.lam, dim / 3.0 * y.mu),
            )
        ),
        list,
        np.array,
    )


@curry
def _check(n_phases, n_phases_other, x_data):
    """Various sanity checks on the data, Elastic modulus and Poissons
    ratio.

    Args:
      n_phases: number of phases in the elastic modulus
      n_phases_other: number of phases in the Poissons ratio
      x_data: the microstructures

    Poissons ratio and elasticy arrays must have equal length


    """
    if n_phases != n_phases_other:
        raise RuntimeError("elastic_modulus and poissons_ratio must be the same length")
    if not issubclass(x_data.dtype.type, np.integer):
        raise TypeError("X must be an integer array")
    if np.max(x_data) >= n_phases or np.min(x_data) < 0:
        raise RuntimeError("X must be between 0 and {N}.".format(N=n_phases - 1))
    if not 3 <= len(x_data.shape) <= 4:
        raise RuntimeError("the shape of x_data is incorrect")


def get_fields(shape, delta_x):
    """Get the fields for the displacement and test function

    Args:
      shape: the shape of the domain
      delta_x: the mesh spacing

    Returns:
      tuple of field variables
    """
    return pipe(
        np.array(shape),
        lambda x: gen_block_mesh(
            x * delta_x, x + 1, np.zeros_like(shape), verbose=False
        ),
        lambda x: Domain("domain", x),
        lambda x: x.create_region("region_all", "all"),
        lambda x: Field.from_args("fu", np.float64, "vector", x, approx_order=2),
        lambda x: (
            FieldVariable("u", "unknown", x),
            FieldVariable("v", "test", x, primary_var_name="u"),
        ),
    )


def _get_material(property_array, domain, delta_x):
    """
    Creates an SfePy material from the material property fields for the
    quadrature points.

    Args:
      property_array: array of the properties with shape (n_x, n_y, n_z, 2)
      domain: the Sfepy domain
      delta_x: the grid spacing

    Returns:
      a SfePy material

    """
    reshape = lambda x: np.ascontiguousarray(x.reshape((x.shape[0], 1, 1)))

    def _material_func_(_, coors, mode=None, **__):
        if mode == "qp":
            return pipe(
                np.empty_like(coors, dtype=int),
                lambda x: np.floor(
                    (coors - domain.get_mesh_bounding_box()[0][None]) / delta_x,
                    x,
                    casting="unsafe",
                ),
                lambda x: x.swapaxes(0, 1),
                tuple,
                lambda x: property_array[x],
                lambda x: dict(
                    lam=reshape(x[..., 0]),
                    mu=reshape(x[..., 1]),
                    D=stiffness_from_lame(
                        domain.get_mesh_bounding_box().shape[1],
                        lam=x[..., 0],
                        mu=x[..., 1],
                    ),
                ),
            )
        return None

    return Material("m", function=Function("material_func", _material_func_))


@curry
def get_term(property_array, delta_x, fields):
    """Get the term

    Args:
      property_array: the spatial array of property values
      delta_x: the grid spacing
      fields: the Sfepy u, v fields

    Returns:
      a new term
    """
    return Term.new(
        "dw_lin_elastic_iso(m.lam, m.mu, v, u)",
        Integral("i", order=4),
        fields[0].field.region,
        m=_get_material(property_array, fields[0].field.region.domain, delta_x),
        u=fields[0],
        v=fields[1],
    )


def subdomain_func(x_points=(), y_points=(), z_points=(), max_x=None):
    """
    Creates a function to mask subdomains in Sfepy.

    Args:
      x_points: tuple of lines or points to be masked in the x-plane
      y_points: tuple of lines or points to be masked in the y-plane
      z_points: tuple of lines or points to be masked in the z-plane

    Returns:
      array of masked location indices

    """

    eps = lambda x: 1e-3 * (x[1, -1] - x[0, -1])

    def np_or(seq):
        if seq:
            return fold((lambda x, y: x | y), seq)
        return True

    @curry
    def flag_it(points, coords, index):
        close = lambda x: (coords[:, index] < (x + eps(coords))) & (
            coords[:, index] > (x - eps(coords))
        )
        return pipe(points, map_(close), list, np_or, lambda x: (len(points) == 0) | x)

    def _func(coords, domain=None):  # pylint: disable=unused-argument
        return pipe(
            (x_points, y_points, z_points),
            enumerate,
            map_(lambda x: flag_it(x[1], coords, x[0])),
            list,
            curry(fold)(lambda x, y: x & y),
            lambda x: (x & (coords[:, 0] < (max_x - eps(coords))))
            if max_x is not None
            else x,
            np.where,
            first,
        )

    return _func


@curry
def get_region_func(max_x, dim_string, domain, name_minmax):
    """Generate the Sfepy region

    Args:
      max_x: max value in the x direction
      dim_string: either "x", "y", or "z"
      domain: the Sfepy domain
      name_minmax: tuple of name and min / max

    Returns:
      the Sfepy region
    """
    return pipe(
        subdomain_func(max_x=max_x, **{dim_string + "_points": (name_minmax[1],)}),
        lambda x: Function(dim_string + name_minmax[0], x),
        lambda x: domain.create_region(
            "region_{0}_{1}".format(dim_string, name_minmax[0]),
            "vertices by {0}".format(x.name),
            "facet",
            functions=Functions([x]),
        ),
    )


def get_bc(max_x_func, domain, dim, bc_dict_func):
    """Get the periodic boundary condition

    Args:
      max_x_func: function for finding the maximum value of x
      domain: the Sfepy domain
      dim: the x, y or z direction
      bc_dict_func: function to generate the bc dict

    Returns:
      the boundary condition and the sfepy function
    """
    dim_dict = lambda x: [
        ("x", per.match_x_plane),
        ("y", per.match_y_plane),
        ("z", per.match_z_plane),
    ][dim][x]
    return pipe(
        domain.get_mesh_bounding_box(),
        lambda x: PeriodicBC(
            "periodic_{0}".format(dim_dict(0)),
            list(
                map_(
                    get_region_func(max_x_func(x), dim_dict(0), domain),
                    zip(("plus", "minus"), x[:, dim][::-1]),
                )
            ),
            bc_dict_func(x),
            match="match_{0}_plane".format(dim_dict(0)),
        ),
        lambda x: (x, Function("match_{0}_plane".format(dim_dict(0)), dim_dict(1))),
    )


@curry
def get_periodic_bc_yz(domain, dim):
    """Get the periodic boundary conditoin in the YZ directions

    Args:
      domain: the Sfepy domain
      dim: the x, y or z directions

    Returns:
      the boundary condition and sfepy function
    """

    return get_bc(
        lambda _: None,
        domain,
        dim,
        lambda x: merge({"u.1": "u.1"}, {"u.2": "u.2"} if x.shape[1] == 3 else dict()),
    )


@curry
def get_periodic_bc_x(domain, dim):
    """Get a periodic buondary condition in the X direction

    Args:
      domain: the Sfepy domain
      dim: the x, y or z directions

    Returns:
      the boundary condition and sfepy function
    """
    return get_bc(
        lambda x: domain.get_mesh_bounding_box()[:, 0][1],
        domain,
        dim,
        lambda x: {"u.0": "u.0"},
    )


def get_periodic_bcs(domain):
    """Get the periodic boundary conditions

    Args:
      domain: the Sfepy domain

    Returns:
      the boundary conditions and sfepy functions
    """
    zipped = lambda x, f: pipe(
        range(x, domain.get_mesh_bounding_box().shape[1]),
        map_(f(domain)),
        lambda x: zip(*x),
        list,
    )

    return pipe(
        (zipped(0, get_periodic_bc_yz), zipped(1, get_periodic_bc_x)),
        lambda x: (Conditions(x[0][0] + x[1][0]), Functions(x[0][1] + x[1][1])),
    )


def get_shift_or_fixed_bcs(domain, points_dict_f, name, x_points_f):
    """Generic function for generating fixed or shift boundary conditions

    Args:
      domain: the sfepy domain
      points_dict_f: function to return displacements
      name: the unique of the boundary condition
      x_points_f: function to return x_points for subdomain_func

    Returns:
      a Sfepy EssentialBC
    """

    def func(min_xyz, max_xyz):
        return pipe(
            dict(z_points=(max_xyz[2], min_xyz[2])) if len(min_xyz) == 3 else dict(),
            lambda x: subdomain_func(
                x_points=x_points_f(min_xyz, max_xyz),
                y_points=(max_xyz[1], min_xyz[1]),
                **x,
            ),
            lambda x: Function(f"{name}_x_points", x),
            lambda x: domain.create_region(
                f"region_{name}_points",
                f"vertices by {name}_x_points",
                "vertex",
                functions=Functions([x]),
            ),
            lambda x: EssentialBC(
                f"{name}_points_BC", x, points_dict_f(min_xyz, max_xyz)
            ),
        )

    return func(domain.get_mesh_bounding_box()[0], domain.get_mesh_bounding_box()[1])


def get_displacement_bcs(domain, macro_strain):
    """Get the shift and fixed BCs.

    The shift BC has the the right top and bottom points in x, y and z
    fixed or displaced.

    The fixed BC has the left top and bottom points in x, y and z
    fixed.

    Args:
      domain: an Sfepy domain
      macro_strain: the macro strain

    Returns:
      the Sfepy boundary conditions

    """
    return Conditions(
        [
            get_shift_or_fixed_bcs(
                domain,
                lambda min_, max_: {"u.0": macro_strain * (max_[0] - min_[0])},
                "shift",
                lambda min_, max_: (max_[0],),
            ),
            get_shift_or_fixed_bcs(
                domain,
                lambda min_, max_: merge(
                    {"u.0": 0.0, "u.1": 0.0}, {"u.2": 0.0} if len(min_) == 3 else dict()
                ),
                "fix",
                lambda min_, max_: (min_[0],),
            ),
        ]
    )


def get_linear_combination_bcs(domain, macro_strain):
    """
    The right nodes are periodic with the left nodes but also displaced.

    Args:
      domain: the Sfepy domain
      macro_strain: the macro strain

    Returns:
      linear combination boundary conditions

    """

    def func(min_xyz, max_xyz):
        def shift_(_, coors, __):
            return np.ones_like(coors[:, 0]) * macro_strain * (max_xyz[0] - min_xyz[0])

        return pipe(
            [("plus", max_xyz[0]), ("minus", min_xyz[0])],
            map_(get_region_func(None, "x", domain)),
            list,
            lambda x: LinearCombinationBC(
                "lcbc",
                x,
                {"u.0": "u.0"},
                Function("match_x_plane", per.match_x_plane),
                "shifted_periodic",
                arguments=(Function("shift", shift_),),
            ),
            lambda x: Conditions([x]),
        )

    return func(domain.get_mesh_bounding_box()[0], domain.get_mesh_bounding_box()[1])


def get_nls(evaluator):
    """Get the non-linear solver

    Args:
      evaluator: the problem evaluator

    Returns:
      the non-linear solver
    """
    return Newton(
        {},
        lin_solver=ScipyDirect({}),
        fun=evaluator.eval_residual,
        fun_grad=evaluator.eval_tangent_matrix,
    )


def get_problem(fields, property_array, delta_x, macro_strain):
    """Build the Sfepy problem

    Args:
      fields: the Sfepy fields
      property_array: lame and mu params for every mesh element
      delta_x: the grid spacing
      macro_strain: the macro strain

    Returns:
      the Sfepy problem
    """

    def func(epbcs, functions, ebcs, lcbcs):
        return pipe(
            fields,
            get_term(property_array, delta_x),
            lambda x: [Equation("balance_of_forces", x)],
            Equations,
            lambda x: Problem("elasticity", equations=x, functions=functions),
            do(lambda x: x.time_update(ebcs=ebcs, epbcs=epbcs, lcbcs=lcbcs)),
            do(lambda x: x.set_solver(get_nls(x.get_evaluator()))),
        )

    return pipe(
        fields[0].field.region.domain,
        lambda x: func(
            *get_periodic_bcs(x),
            get_displacement_bcs(x, macro_strain),
            get_linear_combination_bcs(x, macro_strain),
        ),
    )


def get_stress_strain_data(problem, shape, name):
    """Get the stress or strain data

    Args:
      problem: the sovled Sfepy problem
      shape: the shape of the domain
      name: either 'stress' or 'strain' string

    Returns:
      reshaped stress or strain data
    """
    return pipe(
        dict(
            dims=problem.domain.get_mesh_bounding_box().shape[1],
            args=dict(stress="m.D, u", strain="u")[name],
            name=name,
        ),
        lambda x: "ev_cauchy_{name}.{dims}.region_all({args})".format(**x),
        lambda x: problem.evaluate(x, mode="el_avg", copy_materials=False),
        np.squeeze,
        lambda x: np.reshape(x, (shape + x.shape[-1:])),
    )


def get_displacement(vec, shape):
    """Extract the displacement data

    Args:
      vec: the output from problem.solve()
      shape: the grid shape

    Returns:
      displacement field
    """
    return pipe(
        vec.create_output_dict()["u"].data,
        lambda x: np.reshape(x, (tuple(x + 1 for x in shape) + x.shape[-1:])),
    )


@curry
def get_data(shape, problem, vec):
    """Get the data fields given a solved Sfepy problem

    Args:
      shape: the grid shape
      problem: the solved Sfepy problem
      vec: the output from problem.solve()

    Returns
      the reshaped displacement field
    """

    return (
        get_stress_strain_data(problem, shape, "strain"),
        get_displacement(vec, shape),
        get_stress_strain_data(problem, shape, "stress"),
    )

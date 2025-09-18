from collections.abc import Sequence
from textwrap import dedent

import openmc
from openmc.model.surface_composite import OrthogonalBox, \
    RectangularParallelepiped, RightCircularCylinder, ConicalFrustum
from openmc_mcnp_adapter import mcnp_str_to_model, get_openmc_surfaces
import pytest
from pytest import approx, mark


def convert_surface(mnemonic: str, params: Sequence[float]) -> openmc.Surface:
    """Return an OpenMC surface converted from an MCNP surface

    Parameters
    ----------
    mnemonic
        MCNP surface mnemonic (e.g., "CX")
    params
        Parameters for the surface

    Returns
    -------
    Converted surface

    """
    surface = {'id': 1, 'mnemonic': mnemonic, 'coefficients': params, 'reflective': False}
    surfaces = get_openmc_surfaces([surface], {})
    return surfaces[1]


@mark.parametrize(
    "mnemonic, params, expected_type, attrs",
    [
        ("p",  (4.0,  7.0, -2.5, 1.0), openmc.Plane, ("a", "b", "c", "d")),
        ("px", (5.5,), openmc.XPlane, ("x0",)),
        ("py", (-3.0,), openmc.YPlane, ("y0",)),
        ("pz", ( 1.2,), openmc.ZPlane, ("z0",)),
    ],
)
def test_planes(mnemonic, params, expected_type, attrs):
    surf = convert_surface(mnemonic, params)
    assert isinstance(surf, expected_type)
    for attr, value in zip(attrs, params):
        assert getattr(surf, attr) == approx(value)


def test_plane_from_points():
    # Points defining plane y = x - 1
    coeffs = (1.0, 0.0, 0.0,
              2.0, 1.0, 0.0,
              1.0, 0.0, 1.0)
    surf = convert_surface("p", coeffs)
    assert isinstance(surf, openmc.Plane)
    assert surf.a == approx(1.0)
    assert surf.b == approx(-1.0)
    assert surf.c == approx(0.0)
    assert surf.d == approx(1.0)


@mark.parametrize(
    "mnemonic, params",
    [
        ("so", (2.4,)),
        ("s", (1.0, -2.0, 3.0, 4.5)),
        ("sph", (1.0, -2.0, 3.0, 4.5)),
    ],
)
def test_spheres(mnemonic, params):
    surf = convert_surface(mnemonic, params)

    assert isinstance(surf, openmc.Sphere)

    if mnemonic == "so":
        r = params[0]
        assert surf.x0 == 0.0
        assert surf.y0 == 0.0
        assert surf.z0 == 0.0
        assert surf.r  == approx(r)
    else:
        x0, y0, z0, r = params
        assert surf.x0 == approx(x0)
        assert surf.y0 == approx(y0)
        assert surf.z0 == approx(z0)
        assert surf.r  == approx(r)


@mark.parametrize(
    "mnemonic, params, expected_type, center_attrs",
    [
        ("c/x", ( 1.5, -4.0, 2.2), openmc.XCylinder, ("y0", "z0")),
        ("c/y", (-3.3,  0.0, 1.0), openmc.YCylinder, ("x0", "z0")),
        ("c/z", ( 2.0,  5.5, 0.8), openmc.ZCylinder, ("x0", "y0")),
        ("cx", (2.2,), openmc.XCylinder, ("y0", "z0")),
        ("cy", (1.0,), openmc.YCylinder, ("x0", "z0")),
        ("cz", (0.8,), openmc.ZCylinder, ("x0", "y0")),
    ],
)
def test_cylinders(mnemonic, params, expected_type, center_attrs):
    surf = convert_surface(mnemonic, params)

    assert isinstance(surf, expected_type)

    # center coordinates
    if mnemonic in ("cx", "cy", "cz"):
        assert surf.x0 == 0.0
        assert surf.y0 == 0.0
        assert surf.z0 == 0.0
    else:
        for attr, val in zip(center_attrs, params[:-1]):
            assert getattr(surf, attr) == approx(val)

    # radius (last parameter)
    assert surf.r == approx(params[-1])


@mark.parametrize(
    "mnemonic, params, expected_type, attrs",
    [
        # parallel to axis
        ("k/x",  (0.0, -1.5, 2.0, 1.3), openmc.XCone, ("x0", "y0", "z0", "r2")),
        ("k/y",  (3.0,  1.0, 0.0, 0.8), openmc.YCone, ("x0", "y0", "z0", "r2")),
        ("k/z",  (-2.0, 0.5, 4.1, 2.2), openmc.ZCone, ("x0", "y0", "z0", "r2")),

        # on axis
        ("kx", (-4.1, 2.5),  openmc.XCone, ("x0", "r2")),
        ("ky", (3.2, 1.1),  openmc.YCone, ("y0", "r2")),
        ("kz", (-0.4, 0.7),  openmc.ZCone, ("z0", "r2")),
    ],
)
def test_cones_two_sided(mnemonic, params, expected_type, attrs):
    surf = convert_surface(mnemonic, params)
    assert isinstance(surf, expected_type)
    for attr, val in zip(attrs, params):
        assert getattr(surf, attr) == approx(val)

    if "mnemonic" == "kx":
        assert surf.y0 == 0.0
        assert surf.z0 == 0.0
    elif "mnemonic" == "ky":
        assert surf.x0 == 0.0
        assert surf.z0 == 0.0
    elif "mnemonic" == "kz":
        assert surf.x0 == 0.0
        assert surf.y0 == 0.0


def test_ellipsoid_sq():
    coeffs = (1.0, 2.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.2, 0.0, 0.0)
    surf = convert_surface("sq", coeffs)
    a, b, c, d, e, f, g, x, y, z = coeffs
    assert isinstance(surf, openmc.Quadric)
    assert surf.a == approx(a)
    assert surf.b == approx(b)
    assert surf.c == approx(c)
    assert surf.d == surf.e == surf.f == 0.0
    assert surf.g == approx(2*(d - a*x))
    assert surf.h == approx(2*(e - b*y))
    assert surf.j == approx(2*(f - c*z))
    assert surf.k == approx(a*x*x + b*y*y + c*z*z + 2*(d*x + e*y + f*z) + g)


def test_general_quadric_gq():
    coeffs = (1.0, -0.5, 0.8, 0.0, 1.2, -0.3, 2.0, 0.0, -1.1, 0.6)
    surf = convert_surface("gq", coeffs)

    assert isinstance(surf, openmc.Quadric)

    names = ("a", "b", "c", "d", "e", "f", "g", "h", "j", "k")
    for name, val in zip(names, coeffs):
        assert getattr(surf, name) == approx(val)


@mark.parametrize(
    "mnemonic, expected_type",
    [
        ("tx", openmc.XTorus),
        ("ty", openmc.YTorus),
        ("tz", openmc.ZTorus),
    ],
)
def test_torus(mnemonic, expected_type):
    coeffs = (-2.0,  3.5, 5.0, 1.0, 0.3, 0.2)
    surf = convert_surface(mnemonic, coeffs)
    assert isinstance(surf, expected_type)

    names = ("x0", "y0", "z0", "a", "b", "c")
    for name, val in zip(names, coeffs):
        assert getattr(surf, name) == approx(val)


@mark.parametrize(
    "mnemonic, params, expected_type, attr, value",
    [
        ("x", (-4.0, 3.0), openmc.XPlane, "x0", -4.0),
        ("x", (1.0, 2.0, 1.0, 3.0), openmc.XPlane, "x0", 1.0),
        ("x", (0.0, 1.0, 5.0, 1.0), openmc.XCylinder, "r", 1.0),
        ("y", (6.0, 2.0), openmc.YPlane, "y0", 6.0),
        ("y", (2.5, 3.0, 2.5, 6.0), openmc.YPlane, "y0", 2.5),
        ("y", (0.0, 0.8, -4.0, 0.8), openmc.YCylinder, "r", 0.8),
        ("z", (0.0, 1.0), openmc.ZPlane, "z0", 0.0),
        ("z", (-3.0, 4.4, -3.0, 7.7), openmc.ZPlane, "z0", -3.0),
        ("z", (0.0, 2.2, 9.0, 2.2), openmc.ZCylinder, "r", 2.2),
    ],
)
def test_axisymmetric_surfaces(mnemonic, params, expected_type, attr, value):
    """Test conversion of X/Y/Z surfaces"""
    surf = convert_surface(mnemonic, params)
    assert isinstance(surf, expected_type)
    assert getattr(surf, attr) == approx(value)


def test_box_macrobody():
    coeffs = (0.0, 0.0, 0.0,
              1.0, 0.0, 0.0,
              0.0, 2.0, 0.0,
              0.0, 0.0, 3.0)
    surf = convert_surface("box", coeffs)
    assert isinstance(surf, OrthogonalBox)
    # Plane position along an axis is d / coefficient
    assert surf.ax1_min.d / surf.ax1_min.a == approx(0.0)
    assert surf.ax1_max.d / surf.ax1_max.a == approx(1.0)
    assert surf.ax2_min.d / surf.ax2_min.b == approx(0.0)
    assert surf.ax2_max.d / surf.ax2_max.b == approx(2.0)
    assert surf.ax3_min.d / surf.ax3_min.c == approx(0.0)
    assert surf.ax3_max.d / surf.ax3_max.c == approx(3.0)


def test_rpp_macrobody():
    coeffs = (-1.0, 2.0, -3.0, 4.0, 0.5, 5.5)
    surf = convert_surface("rpp", coeffs)
    assert isinstance(surf, RectangularParallelepiped)
    assert surf.xmin.d / surf.xmin.a == approx(-1.0)
    assert surf.xmax.d / surf.xmax.a == approx(2.0)
    assert surf.ymin.d / surf.ymin.b == approx(-3.0)
    assert surf.ymax.d / surf.ymax.b == approx(4.0)
    assert surf.zmin.d / surf.zmin.c == approx(0.5)
    assert surf.zmax.d / surf.zmax.c == approx(5.5)


@mark.parametrize(
    "coeffs, expected_bottom_z, expected_top_z, r",
    [
        # Base at (0,0,0), height +5 along z
        ((0.0, 0.0, 0.0, 0.0, 0.0, 5.0, 1.5), 0.0, 5.0, 1.5),
        # Negative height vector should be flipped internally (known failing behavior)
        pytest.param((0.0, 0.0, 0.0, 0.0, 0.0, -5.0, 1.0), 0.0, -5.0, 1.0,
                     marks=pytest.mark.xfail(reason="Negative height handling currently broken", strict=False)),
    ],
)
def test_rcc_macrobody(coeffs, expected_bottom_z, expected_top_z, r):
    surf = convert_surface("rcc", coeffs)
    assert isinstance(surf, RightCircularCylinder)
    assert surf.cyl.r == approx(r)
    assert surf.bottom.d / surf.bottom.c == approx(expected_bottom_z)
    assert surf.top.d / surf.top.c == approx(expected_top_z)


def test_trc_macrobody():
    coeffs = (0.0, 0.0, 0.0, 0.0, 0.0, 10.0, 1.0, 2.0)
    surf = convert_surface("trc", coeffs)
    assert isinstance(surf, ConicalFrustum)
    assert surf.plane_bottom.d / surf.plane_bottom.c == approx(0.0)
    assert surf.plane_top.d / surf.plane_top.c == approx(10.0)
    # Check points near boundary
    assert (0.99, 0., 0.01) in -surf
    assert (1.01, 0., 0.01) in +surf
    assert (1.99, 0., 9.99) in -surf
    assert (2.01, 0., 9.99) in +surf
    assert (0., 0., -0.01) in +surf
    assert (0., 0., 10.01) in +surf


def test_rpp_facets():
    mcnp_str = dedent("""
    title
    1  1 -1.0  -1.1 -1.2
    2  1 -1.0  -1.3 -1.4
    3  1 -1.0  -1.5 -1.6

    1  rpp -1.0 2.0 -3.0 4.0 0.5 5.5

    m1   1001.80c  3.0
    """)
    model = mcnp_str_to_model(mcnp_str)
    cells = model.geometry.get_all_cells()
    assert (0., 0., 0.) in cells[1].region
    assert (-2.0, 0., 0) not in cells[1].region
    assert (2.5, 0., 0.) not in cells[1].region
    assert (0., -1.0, 0.) in cells[2].region
    assert (0., -4.0, 0.) not in cells[2].region
    assert (0., 5.0, 0.) not in cells[2].region
    assert (0., 0., 1.0) in cells[3].region
    assert (0., 0., -2.0) not in cells[3].region
    assert (0., 0., 6.0) not in cells[3].region


# Remaining macrobody / complex surfaces not yet implemented in conversion:
# RHP, HEX, REC, ELL, WED, ARB

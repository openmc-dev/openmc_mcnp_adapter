from collections.abc import Sequence

import openmc
from openmc_mcnp_adapter import mcnp_str_to_model
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
    surf_card = "1  " + mnemonic + " " + " ".join(str(x) for x in params)
    mcnp_str = f"title\n1  1 -1.0  -1\n\n{surf_card}\n\nm1   1001.80c  3.0"
    model = mcnp_str_to_model(mcnp_str)
    return model.geometry.get_all_surfaces()[1]


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


# TODO: SQ


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


# TODO: X, Y, Z, and P

# TODO: General plane defined by three points (P)

# TODO: BOX

# TODO: RPP

# TODO: RCC

# TODO: RHP, HEX

# TODO: REC

# TODO: TRC

# TODO: ELL

# TODO: WED

# TODO: ARB

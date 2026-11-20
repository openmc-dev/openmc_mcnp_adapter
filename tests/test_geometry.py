from textwrap import dedent

import numpy as np
from pytest import mark, approx, param, raises
from openmc_mcnp_adapter import mcnp_str_to_model
from openmc_mcnp_adapter.parse import parse_data


@mark.parametrize("whitespace", ["", " ", "\t"])
def test_cell_complement(whitespace):
    # Cell 2 corresponds to r < 2 intersected with z > 0
    mcnp_str = dedent(f"""
    title
    100  1 1.0  +1 : -2
    2    1 1.0  #{whitespace}100

    1  so 2.0
    2  pz 0.0

    m1   1001.80c  1.0
    """)
    model = mcnp_str_to_model(mcnp_str)
    cell = model.geometry.get_all_cells()[2]

    # Check various points
    assert (0., 0., 0.1) in cell.region
    assert (0., 0., -0.1) not in cell.region
    assert (0., 0., 1.99) in cell.region
    assert (0., 0., 2.01) not in cell.region
    assert (1., 1., 1.) in cell.region
    assert (2., 0., 1.) not in cell.region


def test_likenbut():
    mcnp_str = dedent("""
    title
    1   1 -1.0  -1
    2   LIKE 1 BUT MAT=2 RHO=-2.0 TRCL=(2.0 0.0 0.0)

    1   so 1.0

    m1   1001.80c  1.0
    m2   1002.80c  1.0
    """)
    model = mcnp_str_to_model(mcnp_str)
    cell = model.geometry.get_all_cells()[2]

    # Material should be changed to m2
    mat = cell.fill
    assert 'H2' in mat.get_nuclide_densities()

    # Density should be 2.0 g/cm3
    assert mat.get_mass_density() == approx(2.0)

    # Points should correspond to sphere of r=1 centered at (2, 0, 0)
    assert (2.0, 0.0, 0.0) in cell.region
    assert (0.0, 0.0, 0.0) not in cell.region
    assert (2.0, 0.9, 0.0) in cell.region
    assert (2.0, 1.1, 0.0) not in cell.region


@mark.parametrize(
    "cell_card, surface_cards, points_inside, points_outside",
    [
        (
            "1   1 -1.0  -1 TRCL=(2.0 0.0 0.0)",
            ("1   so 1.0",),
            [
                (2.0, 0.0, 0.0),
                (2.0, 0.9, 0.0),
                (2.0, 0.0, 0.9),
            ],
            [
                (0.9, 0.0, 0.0),
                (2.0, 1.1, 0.0),
                (2.0, 0.0, 1.1),
            ],
        ),
        (
            "1 0 -1 TRCL=(1.0 0.0 0.0 0.0 -1.0 0.0 1.0 0.0 0.0 0.0 0.0 1.0)",
            ("1 rpp -0.5 0.5 -0.25 0.25 -0.1 0.1",),
            [
                (1.0, 0.0, 0.0),
                (1.2, 0.0, 0.0),
                (1.0, 0.4, 0.0),
                (1.0, -0.4, 0.0),
            ],
            [
                (0.7, 0.0, 0.0),
                (1.3, 0.0, 0.0),
                (1.0, 0.6, 0.0),
                (1.0, 0.0, 0.2),
            ],
        ),
        (
            "1 0 -1 *TRCL=(1.0 0.0 0.0 90.0 180.0 90.0 0.0 90.0 90.0 90.0 90.0 0.0)",
            ("1 rpp -0.5 0.5 -0.25 0.25 -0.1 0.1",),
            [
                (1.0, 0.0, 0.0),
                (1.2, 0.0, 0.0),
                (1.0, 0.4, 0.0),
                (1.0, -0.4, 0.0),
            ],
            [
                (0.7, 0.0, 0.0),
                (1.3, 0.0, 0.0),
                (1.0, 0.6, 0.0),
                (1.0, 0.0, 0.2),
            ],
        ),
        (
            "1 0 -1 TRCL=(0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0 -1)",
            ("1 rpp -0.5 0.5 -0.25 0.25 -0.1 0.1",),
            [
                (-0.4, 0.1, 0.0),
                (0.0, 0.0, 0.0),
                (0.49, -0.24, -0.05),
            ],
            [
                (-0.6, 0.0, 0.0),
                (0.0, 0.5, 0.0),
                (0.0, 0.0, 0.2),
            ],
        ),
    ],
)
def test_trcl(cell_card, surface_cards, points_inside, points_outside):
    surface_block = "\n".join(surface_cards)
    mcnp_str = dedent(f"""
    title
    {cell_card}

    {surface_block}

    m1   1001.80c  1.0
    """)
    model = mcnp_str_to_model(mcnp_str)
    cell = model.geometry.get_all_cells()[1]

    for point in points_inside:
        assert point in cell.region

    for point in points_outside:
        assert point not in cell.region


def test_trcl_macrobody():
    mcnp_str = dedent("""
    title
    1 0 -1 trcl=(2.0 0.0 0.0)

    1 rpp -1.0 1.0 -1.0 1.0 -1.0 1.0

    m1     1001.80c  1.0
    """)
    model = mcnp_str_to_model(mcnp_str)
    cell = model.geometry.get_all_cells()[1]
    assert (1.5, 0., 0.) in cell.region
    assert (0., 0., 0.) not in cell.region


@mark.parametrize(
    "keywords",
    [
        "FILL=10 TRCL=(2.0 0.0 0.0)",
        "FILL=10(2.0 0.0 0.0)",
        "FILL=10 TRCL=(2.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0)",
        "FILL=10(2.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0)",
        "FILL=10 *TRCL=(2.0 0.0 0.0 0.0 90.0 90.0 90.0 0.0 90.0 90.0 90.0 0.0)",
        "*FILL=10(2.0 0.0 0.0 0.0 90.0 90.0 90.0 0.0 90.0 90.0 90.0 0.0)",
        "FILL=10 TRCL=(2.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0 1)",
        "FILL=10 TRCL=(-2.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0 -1)",
        "FILL=10(2.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0 1)",
        "FILL=10(-2.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0 -1)",
        "FILL=10 *TRCL=(2.0 0.0 0.0 0.0 90.0 90.0 90.0 0.0 90.0 90.0 90.0 0.0 1)",
        "FILL=10 *TRCL=(-2.0 0.0 0.0 0.0 90.0 90.0 90.0 0.0 90.0 90.0 90.0 0.0 -1)",
        "*FILL=10(2.0 0.0 0.0 0.0 90.0 90.0 90.0 0.0 90.0 90.0 90.0 0.0 1)",
        "*FILL=10(-2.0 0.0 0.0 0.0 90.0 90.0 90.0 0.0 90.0 90.0 90.0 0.0 -1)",
        "FILL=10(1)",
        "FILL=10(2)",
        "FILL=10(3)",
    ]
)
def test_fill_transformation(keywords):
    mcnp_str = dedent(f"""
    title
    1 0 -1 {keywords}
    2 0 -2 U=10
    3 0 +2 U=10

    1 so 10.0
    2 so 1.0

    m1     1001.80c  1.0
    tr1    2.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0 1
    *tr2   -2.0 0.0 0.0 0.0 90.0 90.0 90.0 0.0 90.0 90.0 90.0 0.0 -1
    tr3    2.0 0.0 0.0
    """)
    model = mcnp_str_to_model(mcnp_str)
    geometry = model.geometry
    cells = geometry.get_all_cells()

    # Make sure that the cells in universe 10 were shifted
    assert geometry.find((2.0, 0.0, 0.0))[-1] is cells[2]
    assert geometry.find((4.0, 0.0, 0.0))[-1] is cells[3]
    assert geometry.find((0.0, 0.0, 0.0))[-1] is cells[3]


@mark.parametrize(
    "trcl,tr_card",
    [
        # Pure translation along x
        ("TRCL=4", "tr4 2.0 0.0 0.0"),
        # Translation plus identity rotation, cosine and degree forms
        ("TRCL=4", "tr4 2.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0"),
        ("TRCL=4", "*tr4 2.0 0.0 0.0 0.0 90.0 90.0 90.0 0.0 90.0 90.0 90.0 0.0"),
        # The asterisk on TRCL is irrelevant when referring to a TR card
        ("*TRCL=4", "*tr4 2.0 0.0 0.0 0.0 90.0 90.0 90.0 0.0 90.0 90.0 90.0 0.0"),
        # m = -1 reverses the displacement vector
        ("TRCL=4", "tr4 -2.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0 -1"),
    ]
)
def test_trcl_tr_card(trcl, tr_card):
    """A cell TRCL that refers to a TR card by number."""
    mcnp_str = dedent(f"""
    title
    1 0 -1 {trcl}

    1 rpp -1.0 1.0 -1.0 1.0 -1.0 1.0

    m1   1001.80c  1.0
    {tr_card}
    """)
    model = mcnp_str_to_model(mcnp_str)
    cell = model.geometry.get_all_cells()[1]

    # The box should have been shifted to be centered on x = 2
    assert (2.0, 0.0, 0.0) in cell.region
    assert (2.9, 0.9, 0.9) in cell.region
    assert (0.0, 0.0, 0.0) not in cell.region
    assert (4.0, 0.0, 0.0) not in cell.region


def test_trcl_tr_card_rotation():
    """A TR card referenced by TRCL that applies a real (90 deg) rotation."""
    mcnp_str = dedent("""
    title
    1 0 -1 TRCL=7

    1 rpp -4.0 4.0 -1.0 1.0 -1.0 1.0

    m1   1001.80c  1.0
    *tr7 0.0 0.0 0.0 90.0 0.0 90.0 180.0 90.0 90.0 90.0 90.0 0.0
    """)
    model = mcnp_str_to_model(mcnp_str)
    cell = model.geometry.get_all_cells()[1]

    # The box, long in x, should now be long in y
    assert (0.0, 3.0, 0.0) in cell.region
    assert (0.0, -3.0, 0.0) in cell.region
    assert (3.0, 0.0, 0.0) not in cell.region
    assert (-3.0, 0.0, 0.0) not in cell.region


@mark.parametrize(
    "short,full",
    [
        # 60 degree rotation about z, as direction cosines
        (
            "tr8 1.0 2.0 3.0  0.5 0.8660254 0.0  -0.8660254 0.5 0.0",
            "tr8 1.0 2.0 3.0  0.5 0.8660254 0.0  -0.8660254 0.5 0.0  0.0 0.0 1.0",
        ),
        # 30 degree rotation about z, in degrees
        (
            "*tr8 1.0 2.0 3.0  30.0 60.0 90.0  120.0 30.0 90.0",
            "*tr8 1.0 2.0 3.0  30.0 60.0 90.0  120.0 30.0 90.0  90.0 90.0 0.0",
        ),
        # A rotation that does not leave the z-axis fixed: the derived vector
        # is the cross product of the two given ones, here (0, -1, 0)
        (
            "tr8 1.0 2.0 3.0  1.0 0.0 0.0  0.0 0.0 1.0",
            "tr8 1.0 2.0 3.0  1.0 0.0 0.0  0.0 0.0 1.0  0.0 -1.0 0.0",
        ),
        # Same, in degrees
        (
            "*tr8 1.0 2.0 3.0  0.0 90.0 90.0  90.0 90.0 0.0",
            "*tr8 1.0 2.0 3.0  0.0 90.0 90.0  90.0 90.0 0.0  90.0 180.0 90.0",
        ),
        # The m = -1 flag sits right after a 6-entry rotation matrix
        (
            "tr8 -1.0 -2.0 -3.0  0.5 0.8660254 0.0  -0.8660254 0.5 0.0  -1",
            "tr8 -1.0 -2.0 -3.0  0.5 0.8660254 0.0  -0.8660254 0.5 0.0  0.0 0.0 1.0  -1",
        ),
    ]
)
def test_tr_card_nine_values(short, full):
    """A 6-entry rotation matrix is completed with the cross product."""
    def transform(tr_card):
        mcnp_str = dedent(f"""
        title
        1 0 -1 TRCL=8

        1 rpp -4.0 4.0 -1.0 1.0 -1.0 1.0

        m1   1001.80c  1.0
        {tr_card}
        """)
        return mcnp_str_to_model(mcnp_str)

    displacement_short, rotation_short = parse_data(short)['tr'][8]
    displacement_full, rotation_full = parse_data(full)['tr'][8]

    assert displacement_short == approx(displacement_full)
    assert rotation_short == approx(rotation_full)

    # And the resulting geometry agrees
    cell_short = transform(short).geometry.get_all_cells()[1]
    cell_full = transform(full).geometry.get_all_cells()[1]
    for point in [(1.0, 2.0, 3.0), (4.0, 5.0, 3.0), (0.0, 0.0, 0.0), (1.0, 2.0, 5.0)]:
        assert (point in cell_short.region) == (point in cell_full.region)


@mark.parametrize(
    "reversed_card,equivalent",
    [
        # 90 degree rotation about z. The main origin sits at (1, 0, 0) in the
        # auxiliary system, which puts the auxiliary origin at (0, -1, 0) in
        # the main one -- negating alone would wrongly give (-1, 0, 0).
        (
            "tr8 1.0 0.0 0.0  0.0 1.0 0.0  -1.0 0.0 0.0  0.0 0.0 1.0  -1",
            "tr8 0.0 -1.0 0.0  0.0 1.0 0.0  -1.0 0.0 0.0  0.0 0.0 1.0",
        ),
        # Same, with a 6-entry rotation matrix
        (
            "tr8 1.0 0.0 0.0  0.0 1.0 0.0  -1.0 0.0 0.0  -1",
            "tr8 0.0 -1.0 0.0  0.0 1.0 0.0  -1.0 0.0 0.0",
        ),
        # Same, in degrees
        (
            "*tr8 1.0 0.0 0.0  90.0 0.0 90.0  180.0 90.0 90.0  90.0 90.0 0.0  -1",
            "*tr8 0.0 -1.0 0.0  90.0 0.0 90.0  180.0 90.0 90.0  90.0 90.0 0.0",
        ),
        # An identity rotation must keep behaving as a plain negation
        (
            "tr8 -2.0 1.0 0.0  1.0 0.0 0.0  0.0 1.0 0.0  0.0 0.0 1.0  -1",
            "tr8 2.0 -1.0 0.0  1.0 0.0 0.0  0.0 1.0 0.0  0.0 0.0 1.0",
        ),
    ]
)
def test_tr_card_reversed_displacement(reversed_card, equivalent):
    """m = -1 has to rotate the displacement, not just negate it."""
    displacement, rotation = parse_data(reversed_card)['tr'][8]
    displacement_ref, rotation_ref = parse_data(equivalent)['tr'][8]

    assert rotation == approx(rotation_ref)
    assert displacement == approx(displacement_ref)

    # The main origin must land back at the vector given on the m = -1 card
    o_prime = np.array([float(x) for x in reversed_card.split()[1:4]])
    assert rotation @ o_prime + displacement == approx(np.zeros(3))


@mark.parametrize(
    "tr_card",
    [
        # MCNP's implicit-exponent notation, here on the rotation entries
        "tr8 1.0 2.0 3.0  1.0+0 0.0 0.0  0.0 1.0+0 0.0  0.0 0.0 1.0+0",
        "tr8 1.0 2.0 3.0  1.0+0 0.0 0.0  0.0 1.0+0 0.0",
        # and on the displacement
        "tr8 1.0+0 2.0 3.0  1.0 0.0 0.0  0.0 1.0 0.0",
    ]
)
def test_tr_card_implicit_exponent(tr_card):
    """TR card entries may use MCNP's exponent notation without an 'e'."""
    displacement, rotation = parse_data(tr_card)['tr'][8]
    assert displacement == approx([1.0, 2.0, 3.0])
    assert rotation == approx(np.identity(3))


def test_tr_card_six_entry_rotation_is_proper():
    """A 6-entry rotation matrix must not be completed into a singular one."""
    # The two given vectors are x and z; the derived one has to be -y rather
    # than a repeat of z, which would make the matrix singular
    _, rotation = parse_data("tr4 0.0 0.0 0.0  1.0 0.0 0.0  0.0 0.0 1.0")['tr'][4]
    assert np.linalg.det(rotation) == approx(1.0)
    assert rotation.T @ rotation == approx(np.identity(3))


@mark.parametrize(
    "cell_card,extra_cards,match",
    [
        # Composing with a transformation already on the fill
        ("1 0 -1 FILL=2 (1.0 0.0 0.0) TRCL=4", "", "existing FILL"),
        ("1 0 -1 FILL=2 (1.0 0.0 0.0) TRCL=(2.0 0.0 0.0)", "", "existing FILL"),
        # '*FILL' is not picked up by the fill handling
        ("1 0 -1 *FILL=2 (1.0 0.0 0.0) TRCL=4", "", "[*]FILL"),
        ("1 0 -1 *FILL=2 (1.0 0.0 0.0) TRCL=(2.0 0.0 0.0)", "", "[*]FILL"),
    ]
)
def test_trcl_unsupported_fill_combinations(cell_card, extra_cards, match):
    """TRCL combinations that would silently misplace a fill must raise."""
    mcnp_str = dedent(f"""
    title
    {cell_card}
    2 0 -2 u=2

    1 so 10.0
    2 so 1.0

    m1   1001.80c  1.0
    tr4  2.0 0.0 0.0
    {extra_cards}
    """)
    with raises(NotImplementedError, match=match):
        mcnp_str_to_model(mcnp_str)


@mark.parametrize("trcl", ["TRCL=4", "TRCL=(2.0 0.0 0.0)"])
def test_trcl_on_lattice_unsupported(trcl):
    """TRCL on a lattice cell must raise rather than mangle the array spec."""
    mcnp_str = dedent(f"""
    title
    1 0 -1 u=1 lat=1 fill=2 {trcl}
    2 0 -2 u=2
    3 0 -3 fill=1

    1 rpp -1.0 1.0 -1.0 1.0 -1.0 1.0
    2 so 0.5
    3 so 20.0

    m1   1001.80c  1.0
    tr4  2.0 0.0 0.0
    """)
    with raises(NotImplementedError, match="lattice"):
        mcnp_str_to_model(mcnp_str)


#: A 90 degree rotation about z, which makes m = -1 and the 6-entry shorthand
#: behave differently from a plain negation and a fixed z-axis respectively
_ROT9 = "0.0 1.0 0.0  -1.0 0.0 0.0  0.0 0.0 1.0"
_ROT6 = "0.0 1.0 0.0  -1.0 0.0 0.0"


@mark.parametrize(
    "inline,tr_card",
    [
        # m = -1 has to be applied identically either way
        (f"(1.0 0.0 0.0 {_ROT9} -1)", f"tr4 1.0 0.0 0.0 {_ROT9} -1"),
        (f"(1.0 0.0 0.0 {_ROT6} -1)", f"tr4 1.0 0.0 0.0 {_ROT6} -1"),
        # as does the 6-entry rotation shorthand
        (f"(1.0 2.0 3.0 {_ROT6})", f"tr4 1.0 2.0 3.0 {_ROT6}"),
        # and MCNP's implicit-exponent notation
        ("(1.0+0 2.0 3.0)", "tr4 1.0+0 2.0 3.0"),
    ]
)
def test_inline_and_tr_card_transforms_agree(inline, tr_card):
    """An inline TRCL/FILL transform must match the same one on a TR card."""
    def build(trcl, extra):
        return mcnp_str_to_model(dedent(f"""
        title
        1 0 -1 {trcl}
        2 0 -2 u=10

        1 rpp -4.0 4.0 -1.0 1.0 -1.0 1.0
        2 so 0.5

        m1   1001.80c  1.0
        {extra}
        """))

    # The cell region itself
    region_inline = build(f"TRCL={inline}", "").geometry.get_all_cells()[1].region
    region_tr = build("TRCL=4", tr_card).geometry.get_all_cells()[1].region
    for a, b in zip(region_inline.bounding_box, region_tr.bounding_box):
        assert a == approx(b)

    # And the transformation applied to a fill
    cell_inline = build(f"FILL=10 {inline}", "").geometry.get_all_cells()[1]
    cell_tr = build("FILL=10 (4)", tr_card).geometry.get_all_cells()[1]
    assert cell_inline.translation == approx(cell_tr.translation)
    if cell_tr.rotation is None:
        assert cell_inline.rotation is None
    else:
        assert cell_inline.rotation == approx(cell_tr.rotation)


def test_fill_tr_card_undefined():
    """A fill referring to a missing TR card gives a clear error."""
    mcnp_str = dedent("""
    title
    1 0 -1 FILL=10 (9)
    2 0 -2 u=10

    1 so 10.0
    2 so 1.0

    m1   1001.80c  1.0
    """)
    with raises(ValueError, match="TR9"):
        mcnp_str_to_model(mcnp_str)


def test_trcl_tr_card_undefined():
    """Referring to a TR card that does not exist gives a clear error."""
    mcnp_str = dedent("""
    title
    1 0 -1 TRCL=9

    1 so 1.0

    m1   1001.80c  1.0
    """)
    with raises(ValueError, match="TR9"):
        mcnp_str_to_model(mcnp_str)


def test_cell_volume():
    mcnp_str = dedent("""
    title
    1 0 -1  VOL=5.0

    1 so 1.0

    m1   1001.80c  1.0
    """)
    model = mcnp_str_to_model(mcnp_str)
    cell = model.geometry.get_all_cells()[1]
    assert cell.volume == 5.0

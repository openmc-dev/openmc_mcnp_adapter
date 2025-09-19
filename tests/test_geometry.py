from textwrap import dedent

from pytest import mark, approx, param
from openmc_mcnp_adapter import mcnp_str_to_model


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
        param(
            "1 0 -1 TRCL=(0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0 -1)",
            ("1 rpp -0.5 0.5 -0.25 0.25 -0.1 0.1",),
            [],
            [],
            marks=mark.xfail(reason="13-parameter TRCL not yet supported"),
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


def test_trcl_fill():
    mcnp_str = dedent("""
    title
    1 0 -1 FILL=10 TRCL=(2.0 0.0 0.0)
    2 0 -2 U=10
    3 0 +2 U=10

    1 so 10.0
    2 so 1.0

    m1   1001.80c  1.0
    """)
    model = mcnp_str_to_model(mcnp_str)
    geometry = model.geometry
    cells = geometry.get_all_cells()

    # Make sure that the cells in universe 10 were shifted
    assert geometry.find((2.0, 0.0, 0.0))[-1] is cells[2]
    assert geometry.find((4.0, 0.0, 0.0))[-1] is cells[3]
    assert geometry.find((0.0, 0.0, 0.0))[-1] is cells[3]

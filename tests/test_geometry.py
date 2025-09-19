from textwrap import dedent

from pytest import mark, approx
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

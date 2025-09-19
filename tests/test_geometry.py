from textwrap import dedent

from pytest import mark
from openmc_mcnp_adapter import mcnp_str_to_model


@mark.parametrize("whitespace", ["", " ", "\t"])
def test_cell_complement(whitespace):
    # Cell 2 corresponds to r < 2 intersected with z > 0
    mcnp_str = dedent(f"""
    title
    100  1.0  +1 : -2
    2    1.0  #{whitespace}100

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

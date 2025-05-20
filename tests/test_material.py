import textwrap

from openmc_mcnp_adapter import mcnp_str_to_model
from pytest import approx


def test_material_clones():
    mcnp_model = textwrap.dedent("""
        title
        1   1 -1.0    1 -2
        2   1 -2.0    1  2
        3   1 -2.0   -1 -2
        4   1 -1.0   -1  2

        1   px 0.0
        2   py 1.0

        m1   1001.80c  3.0
    """)

    model = mcnp_str_to_model(mcnp_model)
    cells = model.geometry.get_all_cells()
    assert cells[1].fill.id == cells[4].fill.id == 1
    assert cells[2].fill.id == cells[3].fill.id != 1
    assert cells[1].fill.get_mass_density() == approx(1.0)
    assert cells[2].fill.get_mass_density() == approx(2.0)

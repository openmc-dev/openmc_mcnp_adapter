from textwrap import dedent

from pytest import approx
from openmc_mcnp_adapter import mcnp_str_to_model


def test_repeat_shortcut():
    mcnp_str = dedent("""
    title
    1   0 -1

    1   gq  1.0 3r 0.0 5r

    m1   1001.80c  3.0
    """)
    model = mcnp_str_to_model(mcnp_str)
    surf = model.geometry.get_all_surfaces()[1]
    print(surf)

    # Make sure A, B, C, and D parameters are 1.0
    for attr in 'abcd':
        assert getattr(surf, attr) == 1.0

    # Make sure E, F, G, H, J, and K parameters are 0.0
    for attr in 'efghj':
        assert getattr(surf, attr) == 0.0


def test_comments():
    mcnp_str = dedent("""
    title
    c This is a comment line
    1   1 -5.0 -1  $ This is an end-of-line comment

    c surface block
    1   so  3.0  $ Another comment

    m1   1001.80c  3.0  $ Material comment
    """)
    model = mcnp_str_to_model(mcnp_str)
    surf = model.geometry.get_all_surfaces()[1]
    cell = model.geometry.get_all_cells()[1]
    mat = model.materials[0]

    # Sanity checks
    assert surf.r == 3.0
    assert surf.x0 == 0.0
    assert cell.id == 1
    assert str(cell.region) == "-1"
    assert 'H1' in mat.get_nuclide_densities()
    assert mat.get_mass_density() == approx(5.0)

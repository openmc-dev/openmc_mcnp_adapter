from textwrap import dedent

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

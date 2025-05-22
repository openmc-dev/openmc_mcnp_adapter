import textwrap

from openmc_mcnp_adapter import mcnp_str_to_model


def test_vacuum_sphere():
    mcnp_model = textwrap.dedent("""
        title
        1   1 -1.0    -1 imp:n=1
        2   0          1 imp:n=0

        1   so 30.0

        m1   1001.80c  3.0
    """)

    model = mcnp_str_to_model(mcnp_model)
    surf = model.geometry.get_all_surfaces()[1]
    assert surf.boundary_type == 'vacuum'


def test_vacuum_cylinder():
    mcnp_model = textwrap.dedent("""
        title
        1   1 -1.0    -1 2 -3  imp:n=1.0
        2   0          1:-2:3  imp:n=0.0

        1   cz  1.0
        2   pz -1.0
        3   pz  1.0

        m1   1001.80c  3.0
    """)

    model = mcnp_str_to_model(mcnp_model)
    for surf in model.geometry.get_all_surfaces().values():
        assert surf.boundary_type == 'vacuum'


def test_vacuum_box():
    mcnp_model = textwrap.dedent("""
        title
        1   1 -1.0    -1  imp:n=1.00
        2   0          1  imp:n=0.00

        1   rpp -1 1 -1 1 -1 1

        m1   1001.80c  3.0
    """)

    model = mcnp_str_to_model(mcnp_model)
    for surf in model.geometry.get_all_surfaces().values():
        assert surf.boundary_type == 'vacuum'

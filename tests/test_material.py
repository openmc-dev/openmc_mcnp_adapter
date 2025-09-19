import textwrap

import openmc
from openmc_mcnp_adapter import mcnp_str_to_model
from pytest import approx


def convert_material(mat_card: str, density: float, thermal_card: str = "", **kwargs) -> openmc.Material:
    mcnp_model = textwrap.dedent(f"""
        title
        1   1 {density}    1

        1   px 0.0

        {mat_card}
        {thermal_card}
    """)
    model = mcnp_str_to_model(mcnp_model, **kwargs)
    return model.materials[0]


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


def test_material_suffixes():
    # H1 with XS suffix; O16 without suffix
    mat_card = "m1   1001.80c 2.0  8016 1.0"
    thermal_card = "mt1  lwtr grph.10t"
    m = convert_material(mat_card, -2.0, thermal_card)
    nd = m.get_nuclide_densities()
    assert set(nd.keys()) == {'H1', 'O16'}
    assert nd['H1'].percent == approx(2.0)
    assert nd['H1'].percent_type == 'ao'
    assert nd['O16'].percent == approx(1.0)
    assert nd['O16'].percent_type == 'ao'

    # Check S(a,b) tables mapped via get_thermal_name
    # Access private field because there is no public accessor
    sab_names = {name for (name, _) in getattr(m, '_sab', [])}
    assert 'c_H_in_H2O' in sab_names
    assert 'c_Graphite' in sab_names


def test_weight_fractions():
    # 6000 -> natural C, 5010/5011 -> B-10/B-11; negative => weight
    mat_card = "m1   6000 -0.12  5010 -0.2  5011 -0.8"
    m = convert_material(mat_card, -1.0)
    nd = m.get_nuclide_densities()

    # Natural carbon represented as C0 in OpenMC
    assert 'C0' in nd and nd['C0'].percent_type == 'wo' and nd['C0'].percent == approx(0.12)
    assert 'B10' in nd and nd['B10'].percent_type == 'wo' and nd['B10'].percent == approx(0.2)
    assert 'B11' in nd and nd['B11'].percent_type == 'wo' and nd['B11'].percent == approx(0.8)


def test_no_expand_elements():
    # With expand_elements=False, natural oxygen should be added as O0 directly
    mat_card = "m1   47000  1.0"
    m = convert_material(mat_card, -1.0, expand_elements=False)
    nd = m.get_nuclide_densities()
    assert 'Ag0' in nd
    assert nd['Ag0'].percent_type == 'ao'
    assert nd['Ag0'].percent == approx(1.0)


def test_mass_density():
    mat_card = "m1   3006.80c 0.5 3007.80c 0.5"
    m = convert_material(mat_card, -3.0)
    assert m.get_mass_density() == approx(3.0)


def test_atom_density():
    mat_card = "m1   3006.80c 0.5 3007.80c 0.5"
    m = convert_material(mat_card, 0.02)
    nuclide_densities = m.get_nuclide_atom_densities()
    assert sum(nuclide_densities.values()) == approx(0.02)
    assert nuclide_densities['Li6'] == approx(0.01)
    assert nuclide_densities['Li7'] == approx(0.01)

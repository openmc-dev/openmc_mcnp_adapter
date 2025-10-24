# SPDX-FileCopyrightText: 2025 UChicago Argonne, LLC and contributors
# SPDX-License-Identifier: MIT

from textwrap import dedent

from openmc_mcnp_adapter import mcnp_str_to_model, parse_cell


def test_multiparticle_importance():
    """Test that multi-particle importance keywords (e.g., imp:n,p) are parsed correctly."""
    mcnp_str = dedent("""
    title
    1   1 -1.0  -1  imp:n,p=1
    2   0        1  imp:n,p=0

    1   so 1.0

    m1   1001.80c  1.0
    mode n p
    """)
    model = mcnp_str_to_model(mcnp_str)

    # Get cells
    cells = model.geometry.get_all_cells()
    cell1 = cells[1]
    cell2 = cells[2]

    # Check that cells exist and have correct IDs
    assert cell1.id == 1
    assert cell2.id == 2

    # Test parse_cell directly
    line1 = "1   1 -1.0  -1  imp:n,p=1"
    cell_dict1 = parse_cell(line1)
    assert 'imp:n,p' in cell_dict1['parameters']
    assert cell_dict1['parameters']['imp:n,p'] == '1'

    line2 = "2   0       1  imp:n,p=0"
    cell_dict2 = parse_cell(line2)
    assert 'imp:n,p' in cell_dict2['parameters']
    assert cell_dict2['parameters']['imp:n,p'] == '0'


def test_multiparticle_importance_three_particles():
    line = "1   1 -1.0  -1  imp:n,|,e=1"
    cell_dict = parse_cell(line)
    assert 'imp:n,|,e' in cell_dict['parameters']
    assert cell_dict['parameters']['imp:n,|,e'] == '1'


def test_multiparticle_importance_with_other_params():
    """Test that multi-particle importance works with other cell parameters."""
    line = "1   1 -1.0  -1  imp:n,p=1 vol=5.0 u=2"
    cell_dict = parse_cell(line)
    assert 'imp:n,p' in cell_dict['parameters']
    assert cell_dict['parameters']['imp:n,p'] == '1'
    assert cell_dict['parameters']['vol'] == '5.0'
    assert cell_dict['parameters']['u'] == '2'


def test_single_particle_importance_still_works():
    """Ensure single-particle importance still works correctly."""
    line = "1   1 -1.0  -1  imp:n=1"
    cell_dict = parse_cell(line)
    assert 'imp:n' in cell_dict['parameters']
    assert cell_dict['parameters']['imp:n'] == '1'

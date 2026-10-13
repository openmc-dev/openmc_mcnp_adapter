from pathlib import Path
import textwrap
import pytest

from openmc_mcnp_adapter import mcnp_str_to_model, mcnp_to_model
from openmc_mcnp_adapter.parse import expand_read_cards


INPUT_DIR = Path(__file__).with_name("inputs")


def test_read_not_found():
    deck = textwrap.dedent("""
        title
        c The next line points to an invalid file
        read file=/badfile.path
    """)
    with pytest.raises(FileNotFoundError):
        mcnp_str_to_model(deck)


def test_dont_read_comments():
    deck = textwrap.dedent("""
        title
        c The next line would point to an invalid file
        c read file=/badfile.path
        1 0 -1
        
        c The next line would also point to an invalid file
        1 so 1.0  $ read file=/badfile.path
        
        c 
        nps 1
    """)
    mcnp_str_to_model(deck)


def test_read_recursive():
    reference = expand_read_cards(INPUT_DIR / "testReadReference.imcnp")
    trial = expand_read_cards(INPUT_DIR / "testRead.imcnp")
    assert trial == reference


def test_recursive_mcnp_to_model():
    mcnp_to_model(INPUT_DIR / "testRead.imcnp")

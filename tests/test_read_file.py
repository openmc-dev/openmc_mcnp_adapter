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


def test_read_recursive():
    reference = expand_read_cards(INPUT_DIR / "testReadReference.imcnp")
    trial = expand_read_cards(INPUT_DIR / "testRead.imcnp")
    assert trial == reference


def test_recursive_mcnp_to_model():
    mcnp_to_model(INPUT_DIR / "testRead.imcnp")

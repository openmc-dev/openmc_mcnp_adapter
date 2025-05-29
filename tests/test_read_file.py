import os
import textwrap
import pytest

from openmc_mcnp_adapter import mcnp_str_to_model
from openmc_mcnp_adapter.parse import read_file

here = os.path.dirname(os.path.abspath(__file__))
input_dir = os.path.join(here, "inputs")


def test_read_not_found():
	deck = textwrap.dedent("""
		title
		c The next line points to an invalid file
		read file=/badfile.path
	""")
	with pytest.raises(FileNotFoundError):
		mcnp_str_to_model(deck)


def test_read_recursive():
	reference = read_file(os.path.join(input_dir, "testReadReference.imcnp"))
	trial = read_file(os.path.join(input_dir, "testRead.imcnp"))
	assert trial == reference

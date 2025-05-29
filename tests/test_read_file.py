import textwrap
import pytest

from openmc_mcnp_adapter import mcnp_str_to_model


def test_read_not_found():
	deck = textwrap.dedent("""
		title
		c The next line points to an invalid file
		read file=/badfile.path
	""")
	with pytest.raises(FileNotFoundError):
		mcnp_str_to_model(deck)

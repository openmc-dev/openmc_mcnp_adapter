# MCNP Conversion Tools for OpenMC

[![License](https://img.shields.io/badge/license-MIT-green)](https://opensource.org/licenses/MIT)
[![codecov](https://codecov.io/github/openmc-dev/openmc_mcnp_adapter/graph/badge.svg?token=KX00MQ57G5)](https://codecov.io/github/openmc-dev/openmc_mcnp_adapter)

This repository provides tools for parsing/converting MCNP models to OpenMC
classes and/or XML files. To install these tools, run:

    python -m pip install git+https://github.com/openmc-dev/openmc_mcnp_adapter.git

This makes the `openmc_mcnp_adapter` Python module and `mcnp_to_openmc` console
script available. To convert an MCNP model, run:

    mcnp_to_openmc mcnp_input

## Disclaimer

There has been no methodical V&V on this converter; use at your own risk!

## Known Limitations

The converter currently only handles geometry and material information; source
definition (SDEF) and tally specifications are ignored.

The converter will try to set surface boundary conditions to match the MCNP
model, but in many cases it doesn't work cleanly. For these cases, you will need
to manually set boundary conditions on the outermost surfaces.

Some geometry features are not currently supported:

- `X`, `Y`, and `Z` surfaces with 3 coordinate pairs
- `RHP`, `REC`, `ELL`, `WED`, and `ARB` macrobodies
- Hexagonal lattices
- One-dimensional lattices
- Two-dimensional lattices with basis other than x-y
- `U`, `LAT`, and `FILL` cards specified in the data card block

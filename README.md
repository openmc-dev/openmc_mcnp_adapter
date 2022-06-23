# MCNP Conversion Tools for OpenMC

[![License](https://img.shields.io/badge/license-MIT-green)](https://opensource.org/licenses/MIT)

This repository provides tools for parsing/converting MCNP models to OpenMC
classes and/or XML files. To install these tools, run:

    python -m pip install git+https://github.com/openmc-dev/openmc_mcnp_adapter.git

This makes the `openmc_mcnp_adapter` Python module and `mcnp_to_openmc` console
script available. To convert an MCNP model, run:

    mcnp_to_openmc mcnp_input

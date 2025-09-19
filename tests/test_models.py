from pathlib import Path
import sys

import pytest
from openmc_mcnp_adapter.openmc_conversion import mcnp_to_openmc


MODELS_DIR = Path(__file__).parent / "models"
MODELS = sorted(MODELS_DIR.glob("*.mcnp"))


@pytest.mark.parametrize("mcnp_model", MODELS, ids=[path.stem for path in MODELS])
def test_mcnp_models_convert(tmp_path, monkeypatch, mcnp_model):
    output_path = tmp_path / "model.xml"
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, "argv", [
        "mcnp_to_openmc",
        str(mcnp_model),
        "-o",
        str(output_path),
    ])

    mcnp_to_openmc()

    assert output_path.exists(), "Expected OpenMC model.xml was not created"

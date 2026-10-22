from odmantic import Model
from molecular_qm_util.rdkit_scripts.rdkit_optimize import RDKitParams, RDKitForceField
from simstack.models.simstack_model import is_simstack_model


def test_rdkit_params_registration_and_defaults():
    assert issubclass(RDKitParams, Model)
    assert is_simstack_model(RDKitParams)
    params = RDKitParams()
    assert params.max_iters == 500
    assert params.forcefield == RDKitForceField.MMFF94
    assert params.mmff_variant == "MMFF94"
    assert "max_iters" in RDKitParams.json_schema()["properties"]


def test_rdkit_params_round_trip():
    params = RDKitParams(max_iters=12, forcefield="uff", mmff_variant="MMFF94s")
    restored = RDKitParams.model_validate_doc(params.model_dump_doc())
    assert restored.max_iters == 12
    assert restored.forcefield == "uff"
    assert restored.mmff_variant == "MMFF94s"

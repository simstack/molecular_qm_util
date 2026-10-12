import copy

import numpy as np
import pytest

from molecular_qm_models import MoleculeList
from molecular_qm_util import (
    RDKitForceFieldParams,
    RDKitOptimizationParams,
    evaluate_molecules_rdkit,
    optimize_molecules_rdkit,
    score_molecules_rdkit,
    smiles_to_molecule,
)


@pytest.mark.parametrize("forcefield", ["mmff", "mmff94s", "uff", "rdkit_uff"])
def test_batch_scoring_preserves_geometry_order_and_metadata(forcefield):
    inputs = [smiles_to_molecule("CCO"), smiles_to_molecule("CCCC")]
    for index, molecule in enumerate(inputs):
        molecule.properties["label"] = {"index": index}
        molecule.properties["energy"] = 9999
    original = copy.deepcopy(inputs)
    scored = evaluate_molecules_rdkit(inputs, forcefield, threads=2)
    optimized = evaluate_molecules_rdkit(inputs, forcefield, optimize=True, max_iters=200, threads=2)
    for index, (source, score, result) in enumerate(zip(inputs, scored, optimized)):
        assert [a.position for a in score.atoms] == [a.position for a in source.atoms]
        assert source.model_dump() == original[index].model_dump()
        assert score.properties["label"] == {"index": index}
        assert result.properties["label"] == {"index": index}
        assert result.properties["energy"] <= score.properties["energy"] + 1e-6
        assert np.isfinite(result.properties["energy"])
        assert result.properties["energy_unit"] == "kcal/mol"
        assert [a.element for a in result.atoms] == [a.element for a in source.atoms]
        assert result.smiles == source.smiles
        score.properties["label"]["index"] = -1
        assert source.properties["label"]["index"] == index


def test_invalid_forcefield_fails_instead_of_falling_back_to_uff():
    with pytest.raises(ValueError, match="Unsupported force field"):
        evaluate_molecules_rdkit([], "typo")


def test_empty_batches():
    assert evaluate_molecules_rdkit([]) == []
    assert evaluate_molecules_rdkit([], optimize=True) == []


def test_missing_parameters_are_not_zero_energy(monkeypatch):
    import importlib
    module = importlib.import_module("molecular_qm_util.rdkit_scripts.rdkit_optimize")
    molecule = smiles_to_molecule("CCO")
    monkeypatch.setattr(module.AllChem, "MMFFGetMoleculeProperties", lambda *a, **kw: None)
    with pytest.raises(ValueError, match="parameters are unavailable"):
        evaluate_molecules_rdkit([molecule])


def test_batch_nodes_score_and_optimize_a_loaded_list():
    molecule = smiles_to_molecule("CCO")
    molecules = MoleculeList()
    molecules.add_molecule(molecule)
    scored = score_molecules_rdkit.__wrapped__(
        molecules, RDKitForceFieldParams(forcefield="uff", threads=1)
    )
    optimized = optimize_molecules_rdkit.__wrapped__(
        molecules, RDKitOptimizationParams(forcefield="uff", threads=1, max_iters=20)
    )
    assert getattr(score_molecules_rdkit, "_is_node", False) is True
    assert getattr(optimize_molecules_rdkit, "_is_node", False) is True
    assert len(scored) == 1
    assert scored[0].properties["energy_unit"] == "kcal/mol"
    assert [atom.position for atom in scored[0].atoms] == [atom.position for atom in molecule.atoms]
    assert optimized[0].properties["energy"] <= scored[0].properties["energy"] + 1e-6


def test_batch_node_rejects_an_unloaded_molecule_list():
    molecules = MoleculeList()
    molecules.elements.append(smiles_to_molecule("CCO").id)
    with pytest.raises(ValueError, match="loaded MoleculeList"):
        score_molecules_rdkit.__wrapped__(molecules, RDKitForceFieldParams())


def test_existing_single_molecule_node_honors_iteration_parameter():
    from molecular_qm_util.rdkit_scripts.rdkit_optimize import RDKitParams, rdkit_optimize

    molecule = smiles_to_molecule("CCO")
    params = RDKitParams()
    params.max_iters = 10
    result = rdkit_optimize.__wrapped__(molecule, params)
    assert np.isfinite(result.properties["energy"])

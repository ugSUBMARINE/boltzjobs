"""Tests for Boltz API input export."""

from __future__ import annotations

import base64

import pytest

from boltzjobs import Job


def test_boltz_api_entity_conversion_with_modifications_and_ligands() -> None:
    job = Job("api_entities")
    protein = job.add_protein_chain("MVTPEGNVSLVDESLLVGVTDED", ids="A", cyclic=True)
    protein.add_modification("MSE", 2)
    dna = job.add_dna_chain("ATCG", ids="D")
    dna.add_modification("5CM", 3)
    rna = job.add_rna_chain("AUCG", ids="R")
    rna.add_modification("PSU", 4)
    job.add_ligand(ccd="ATP", ids="L1")
    job.add_ligand(smiles="CCO", ids="L2")

    api_input = job.to_boltz_api_input()

    assert api_input == {
        "entities": [
            {
                "type": "protein",
                "value": "MVTPEGNVSLVDESLLVGVTDED",
                "chain_ids": ["A"],
                "cyclic": True,
                "modifications": [
                    {"residue_index": 1, "type": "ccd", "value": "MSE"}
                ],
            },
            {
                "type": "dna",
                "value": "ATCG",
                "chain_ids": ["D"],
                "modifications": [
                    {"residue_index": 2, "type": "ccd", "value": "5CM"}
                ],
            },
            {
                "type": "rna",
                "value": "AUCG",
                "chain_ids": ["R"],
                "modifications": [
                    {"residue_index": 3, "type": "ccd", "value": "PSU"}
                ],
            },
            {"type": "ligand_ccd", "value": "ATP", "chain_ids": ["L1"]},
            {"type": "ligand_smiles", "value": "CCO", "chain_ids": ["L2"]},
        ]
    }


def test_boltz_api_msa_conversion_for_empty_local_and_url(tmp_path) -> None:
    a3m_path = tmp_path / "seq.a3m"
    a3m_path.write_text(">A\nMVT\n")

    job = Job("api_msa")
    job.add_protein_chain("MVT", ids="A")
    job.add_protein_chain("MVT", ids="B", msa="empty")
    job.add_protein_chain("MVT", ids="C", msa=str(a3m_path))
    job.add_protein_chain("MVT", ids="D", msa="https://example.org/paired.csv")

    entities = job.to_boltz_api_input()["entities"]

    assert "msa" not in entities[0]
    assert entities[1]["msa"] == {"type": "empty"}
    assert entities[2]["msa"] == {
        "type": "custom",
        "format": "a3m",
        "source": {
            "type": "base64",
            "data": base64.b64encode(b">A\nMVT\n").decode("ascii"),
            "media_type": "text/x-a3m",
        },
    }
    assert entities[3]["msa"] == {
        "type": "custom",
        "format": "csv",
        "source": {"type": "url", "url": "https://example.org/paired.csv"},
    }


def test_boltz_api_constraints_and_bonds_are_zero_indexed() -> None:
    job = Job("api_constraints")
    job.add_protein_chain("MVTPEGNVSLVDESLLVGVTDED", ids="A")
    job.add_ligand(ccd="ATP", ids="L")
    job.add_bond("A", 2, "CA", "L", 1, "C1")
    job.add_contact("A", 3, "L", "C1", max_distance=8.0, force=True)
    pocket = job.add_pocket("L", max_distance=6.0)
    pocket.add_contact_token("A", 4)
    pocket.add_contact_token("A", 8)

    api_input = job.to_boltz_api_input()

    assert api_input["bonds"] == [
        {
            "atom1": {
                "type": "polymer_atom",
                "chain_id": "A",
                "residue_index": 1,
                "atom_name": "CA",
            },
            "atom2": {"type": "ligand_atom", "chain_id": "L", "atom_name": "C1"},
        }
    ]
    assert api_input["constraints"] == [
        {
            "type": "contact",
            "token1": {"type": "polymer_contact", "chain_id": "A", "residue_index": 2},
            "token2": {"type": "ligand_contact", "chain_id": "L", "atom_name": "C1"},
            "max_distance_angstrom": 8.0,
            "force": True,
        },
        {
            "type": "pocket",
            "binder_chain_id": "L",
            "contact_residues": {"A": [3, 7]},
            "max_distance_angstrom": 6.0,
        },
    ]


def test_boltz_api_templates_support_base64_url_defaults_and_force(tmp_path) -> None:
    cif_path = tmp_path / "template.cif"
    cif_path.write_text("data_test\n")

    job = Job("api_templates")
    job.add_protein_chain("MVTPEGNVSLVDESLLVGVTDED", ids="A")
    job.add_template(cif=str(cif_path), chain_id="A")
    job.add_template(
        pdb="https://example.org/template.pdb",
        chain_id=["A"],
        template_id=["B"],
        force=True,
        threshold=2.5,
    )

    templates = job.to_boltz_api_input()["templates"]

    assert templates[0] == {
        "template_structure": {
            "type": "base64",
            "data": base64.b64encode(b"data_test\n").decode("ascii"),
            "media_type": "chemical/x-cif",
        },
        "template_chains": [{"input_chain_id": "A", "template_chain_id": "A"}],
    }
    assert templates[1] == {
        "template_structure": {
            "type": "url",
            "url": "https://example.org/template.pdb",
        },
        "template_chains": [{"input_chain_id": "A", "template_chain_id": "B"}],
        "force_threshold_angstroms": 2.5,
    }


def test_boltz_api_ligand_binding_from_affinity() -> None:
    job = Job("api_ligand_binding")
    job.add_protein_chain("MVTPEGNVSLVDESLLVGVTDED", ids="A")
    job.add_ligand(smiles="CCO", ids="L")
    job.request_affinity("L")

    api_input = job.to_boltz_api_input(num_samples=3)

    assert api_input["binding"] == {
        "type": "ligand_protein_binding",
        "binder_chain_id": "L",
    }
    assert api_input["num_samples"] == 3


def test_boltz_api_protein_binding_from_method_argument() -> None:
    job = Job("api_protein_binding")
    job.add_protein_chain("MVTPEGNVSLVDESLLVGVTDED", ids="A")
    job.add_protein_chain("MVTPEGNVSLVDESLLVGVTDED", ids=["B", "C"])

    api_input = job.to_boltz_api_input(
        protein_binder_chain_ids=["B", "C"],
        recycling_steps=4,
        sampling_steps=100,
        step_scale=1.5,
    )

    assert api_input["binding"] == {
        "type": "protein_protein_binding",
        "binder_chain_ids": ["B", "C"],
    }
    assert api_input["model_options"] == {
        "recycling_steps": 4,
        "sampling_steps": 100,
        "step_scale": 1.5,
    }


@pytest.mark.parametrize(
    ("option_name", "value", "message"),
    [
        ("num_samples", 0, "num_samples"),
        ("num_samples", 11, "num_samples"),
        ("recycling_steps", 0, "recycling_steps"),
        ("sampling_steps", 49, "sampling_steps"),
        ("step_scale", 1.2, "step_scale"),
        ("step_scale", 2.1, "step_scale"),
    ],
)
def test_boltz_api_option_validation(
    option_name: str, value: int | float, message: str
) -> None:
    job = Job("api_option_validation")
    job.add_protein_chain("MVTPEGNVSLVDESLLVGVTDED", ids="A")

    with pytest.raises(ValueError, match=message):
        if option_name == "num_samples":
            job.to_boltz_api_input(num_samples=int(value))
        elif option_name == "recycling_steps":
            job.to_boltz_api_input(recycling_steps=int(value))
        elif option_name == "sampling_steps":
            job.to_boltz_api_input(sampling_steps=int(value))
        elif option_name == "step_scale":
            job.to_boltz_api_input(step_scale=float(value))
        else:
            raise AssertionError(f"Unexpected option: {option_name}")


def test_boltz_api_missing_local_file_raises() -> None:
    job = Job("api_missing_file")
    job.add_protein_chain("MVTPEGNVSLVDESLLVGVTDED", ids="A", msa="missing.a3m")

    with pytest.raises(ValueError, match="Local file does not exist"):
        job.to_boltz_api_input()


def test_boltz_api_unsupported_msa_extension_raises(tmp_path) -> None:
    msa_path = tmp_path / "seq.txt"
    msa_path.write_text("MVT")
    job = Job("api_bad_msa")
    job.add_protein_chain("MVTPEGNVSLVDESLLVGVTDED", ids="A", msa=str(msa_path))

    with pytest.raises(ValueError, match="Unsupported MSA extension"):
        job.to_boltz_api_input()


def test_boltz_api_ligand_binding_rejects_multi_copy_ligand() -> None:
    job = Job("api_multi_ligand")
    job.add_protein_chain("MVTPEGNVSLVDESLLVGVTDED", ids="A")
    job.add_ligand(smiles="CCO", ids=["L1", "L2"])
    job.request_affinity("L1")

    with pytest.raises(ValueError, match="exactly one chain ID"):
        job.to_boltz_api_input()


def test_boltz_api_ligand_binding_rejects_nucleic_acids() -> None:
    job = Job("api_nucleic_binding")
    job.add_protein_chain("MVTPEGNVSLVDESLLVGVTDED", ids="A")
    job.add_dna_chain("ATCG", ids="D")
    job.add_ligand(smiles="CCO", ids="L")
    with pytest.warns(UserWarning, match="Non-protein targets detected"):
        job.request_affinity("L")

    with pytest.raises(ValueError, match="only supports complexes containing proteins and ligands"):
        job.to_boltz_api_input()


def test_boltz_api_rejects_conflicting_binding_modes() -> None:
    job = Job("api_conflicting_binding")
    job.add_protein_chain("MVTPEGNVSLVDESLLVGVTDED", ids="A")
    job.add_protein_chain("MVTPEGNVSLVDESLLVGVTDED", ids="B")
    job.add_ligand(smiles="CCO", ids="L")
    job.request_affinity("L")

    with pytest.raises(ValueError, match="both ligand-protein affinity and protein-protein"):
        job.to_boltz_api_input(protein_binder_chain_ids="B")


def test_boltz_api_rejects_invalid_protein_binder_chain() -> None:
    job = Job("api_invalid_protein_binder")
    job.add_protein_chain("MVTPEGNVSLVDESLLVGVTDED", ids="A")
    job.add_ligand(ccd="ATP", ids="L")

    with pytest.raises(ValueError, match="requires protein binder chains"):
        job.to_boltz_api_input(protein_binder_chain_ids="L")


def test_boltz_api_rejects_smiles_ligand_atom_references() -> None:
    job = Job("api_smiles_atom")
    job.add_protein_chain("MVTPEGNVSLVDESLLVGVTDED", ids="A")
    job.add_ligand(smiles="CCO", ids="L")
    job.add_contact("A", 1, "L", "C1")

    with pytest.raises(ValueError, match="SMILES ligand chain"):
        job.to_boltz_api_input()


def test_boltz_api_rejects_pocket_ligand_atom_contacts() -> None:
    job = Job("api_bad_pocket")
    job.add_protein_chain("MVTPEGNVSLVDESLLVGVTDED", ids="A")
    job.add_ligand(ccd="ATP", ids="L")
    pocket = job.add_pocket("L")
    pocket.add_contact_token("L", "C1")

    with pytest.raises(ValueError, match="Pocket contacts"):
        job.to_boltz_api_input()


def test_boltz_api_template_validation_errors(tmp_path) -> None:
    cif_path = tmp_path / "template.cif"
    cif_path.write_text("data_test\n")

    no_chain = Job("api_template_no_chain")
    no_chain.add_protein_chain("MVTPEGNVSLVDESLLVGVTDED", ids="A")
    no_chain.add_template(cif=str(cif_path))
    with pytest.raises(ValueError, match="require at least one chain_id"):
        no_chain.to_boltz_api_input()

    mismatch = Job("api_template_mismatch")
    mismatch.add_protein_chain("MVTPEGNVSLVDESLLVGVTDED", ids=["A", "B"])
    mismatch.add_template(cif=str(cif_path), chain_id=["A", "B"], template_id=["X"])
    with pytest.raises(ValueError, match="matching chain_id and template_id lengths"):
        mismatch.to_boltz_api_input()

    too_many = Job("api_template_too_many")
    too_many.add_protein_chain("MVTPEGNVSLVDESLLVGVTDED", ids="A")
    for idx in range(5):
        too_many.add_template(cif=str(cif_path), chain_id="A", template_id=f"T{idx}")
    with pytest.raises(ValueError, match="at most 4 templates"):
        too_many.to_boltz_api_input()


def test_boltz_api_empty_job_raises() -> None:
    with pytest.raises(ValueError, match="Empty list of sequences"):
        Job("empty").to_boltz_api_input()

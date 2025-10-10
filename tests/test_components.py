"""Tests for boltzjobs.components module."""

import pytest
from dataclasses import FrozenInstanceError

from boltzjobs.components import (
    ProteinChain,
    DnaChain,
    RnaChain,
    Ligand,
    Bond,
    Contact,
    Pocket,
    Template,
    Affinity,
    SequenceModification,
    Chain,
)
from boltzjobs.utils import FlowStyleList, SingleQuoted


class TestSequenceModification:
    """Test SequenceModification class."""

    @pytest.mark.unit
    def test_sequence_modification_creation(self):
        """Test basic SequenceModification creation."""
        mod = SequenceModification("HIS", 5)
        assert mod.ccd == "HIS"
        assert mod.position == 5

    @pytest.mark.unit
    def test_sequence_modification_str(self):
        """Test SequenceModification string representation."""
        mod = SequenceModification("CYS", 10)
        assert str(mod) == "   Modification: CYS at position 10"

    @pytest.mark.unit
    def test_sequence_modification_from_dict(self):
        """Test SequenceModification.from_dict()."""
        data = {"ccd": "TYR", "position": 15}
        mod = SequenceModification.from_dict(data)
        assert mod.ccd == "TYR"
        assert mod.position == 15

    @pytest.mark.unit
    def test_sequence_modification_from_dict_defaults(self):
        """Test SequenceModification.from_dict() with defaults."""
        mod = SequenceModification.from_dict({})
        assert mod.ccd == ""
        assert mod.position == 1


class TestChain:
    """Test base Chain class."""

    @pytest.mark.unit
    def test_chain_creation(self, sample_protein_sequence):
        """Test basic Chain creation."""
        chain = Chain(["A"], sample_protein_sequence)
        assert chain.ids == ["A"]
        assert chain.sequence == sample_protein_sequence
        assert chain.modifications == []
        assert chain.cyclic is False

    @pytest.mark.unit
    def test_chain_empty_sequence_error(self):
        """Test Chain creation with empty sequence raises error."""
        with pytest.raises(ValueError, match="Sequence cannot be empty"):
            Chain(["A"], "")

    @pytest.mark.unit
    def test_chain_add_modification(self, sample_protein_sequence):
        """Test adding modifications to chain."""
        chain = Chain(["A"], sample_protein_sequence)
        result = chain.add_modification("HIS", 5)

        assert len(chain.modifications) == 1
        assert chain.modifications[0].ccd == "HIS"
        assert chain.modifications[0].position == 5
        assert result is chain  # Should return self for chaining

    @pytest.mark.unit
    def test_chain_set_cyclic(self, sample_protein_sequence):
        """Test setting cyclic property."""
        chain = Chain(["A"], sample_protein_sequence)
        chain.set_cyclic(True)
        assert chain.cyclic is True

    @pytest.mark.unit
    def test_chain_set_cyclic_invalid_type(self, sample_protein_sequence):
        """Test setting cyclic with invalid type."""
        chain = Chain(["A"], sample_protein_sequence)
        with pytest.raises(TypeError, match="Cyclic must be a boolean value"):
            chain.set_cyclic("true")

    @pytest.mark.unit
    def test_chain_str_short_sequence(self):
        """Test Chain string representation with short sequence."""
        chain = Chain(["A"], "MVTP")
        str_repr = str(chain)
        assert "ID: A" in str_repr
        assert "Sequence: MVTP, 4 residues" in str_repr

    @pytest.mark.unit
    def test_chain_str_long_sequence(self):
        """Test Chain string representation with long sequence."""
        long_seq = "A" * 50
        chain = Chain(["A"], long_seq)
        str_repr = str(chain)
        assert "AAAAAAAAAA.....AAAAAAAAAA" in str_repr

    @pytest.mark.unit
    def test_chain_str_multiple_ids(self):
        """Test Chain string representation with multiple IDs."""
        chain = Chain(["A", "B", "C"], "MVTP")
        str_repr = str(chain)
        assert "IDs: A, B, C" in str_repr

    @pytest.mark.unit
    def test_chain_str_with_modifications_and_cyclic(self, sample_protein_sequence):
        """Test Chain string representation with modifications and cyclic."""
        chain = Chain(["A"], sample_protein_sequence, cyclic=True)
        chain.add_modification("HIS", 5)
        str_repr = str(chain)
        assert "Cyclic: True" in str_repr
        assert "Modification: HIS at position 5" in str_repr

    @pytest.mark.unit
    def test_chain_to_dict_basic(self, sample_protein_sequence):
        """Test Chain.to_dict() basic functionality."""
        chain = Chain(["A"], sample_protein_sequence)
        result = chain.to_dict()

        expected = {"id": "A", "sequence": sample_protein_sequence}
        assert result == expected

    @pytest.mark.unit
    def test_chain_to_dict_multiple_ids(self, sample_protein_sequence):
        """Test Chain.to_dict() with multiple IDs."""
        chain = Chain(["A", "B"], sample_protein_sequence)
        result = chain.to_dict()

        assert isinstance(result["id"], FlowStyleList)
        assert result["id"] == ["A", "B"]

    @pytest.mark.unit
    def test_chain_to_dict_with_modifications_and_cyclic(self, sample_protein_sequence):
        """Test Chain.to_dict() with modifications and cyclic."""
        chain = Chain(["A"], sample_protein_sequence, cyclic=True)
        chain.add_modification("HIS", 5)
        chain.add_modification("CYS", 10)

        result = chain.to_dict()

        assert result["cyclic"] is True
        assert len(result["modifications"]) == 2
        assert result["modifications"][0] == {"position": 5, "ccd": "HIS"}
        assert result["modifications"][1] == {"position": 10, "ccd": "CYS"}

    @pytest.mark.unit
    def test_chain_from_dict_basic(self, sample_protein_sequence):
        """Test Chain.from_dict() basic functionality."""
        data = {"id": "A", "sequence": sample_protein_sequence}
        chain = Chain.from_dict(data)

        assert chain.ids == ["A"]
        assert chain.sequence == sample_protein_sequence
        assert chain.cyclic is False
        assert chain.modifications == []

    @pytest.mark.unit
    def test_chain_from_dict_multiple_ids(self, sample_protein_sequence):
        """Test Chain.from_dict() with multiple IDs."""
        data = {"id": ["A", "B"], "sequence": sample_protein_sequence}
        chain = Chain.from_dict(data)

        assert chain.ids == ["A", "B"]

    @pytest.mark.unit
    def test_chain_from_dict_with_modifications_and_cyclic(
        self, sample_protein_sequence
    ):
        """Test Chain.from_dict() with modifications and cyclic."""
        data = {
            "id": "A",
            "sequence": sample_protein_sequence,
            "cyclic": True,
            "modifications": [
                {"position": 5, "ccd": "HIS"},
                {"position": 10, "ccd": "CYS"},
            ],
        }
        chain = Chain.from_dict(data)

        assert chain.cyclic is True
        assert len(chain.modifications) == 2
        assert chain.modifications[0].ccd == "HIS"
        assert chain.modifications[0].position == 5


class TestProteinChain:
    """Test ProteinChain class."""

    @pytest.mark.unit
    def test_protein_chain_creation(self, sample_protein_sequence):
        """Test basic ProteinChain creation."""
        chain = ProteinChain(["A"], sample_protein_sequence)
        assert chain.ids == ["A"]
        assert chain.sequence == sample_protein_sequence
        assert chain.msa is None

    @pytest.mark.unit
    def test_protein_chain_with_msa(self, sample_protein_sequence):
        """Test ProteinChain creation with MSA."""
        chain = ProteinChain(["A"], sample_protein_sequence, msa="test.a3m")
        assert chain.msa == "test.a3m"

    @pytest.mark.unit
    def test_protein_chain_invalid_sequence(self, invalid_protein_sequence):
        """Test ProteinChain with invalid amino acids."""
        with pytest.raises(
            ValueError, match="Protein sequence contains invalid one-letter codes"
        ):
            ProteinChain(["A"], invalid_protein_sequence)

    @pytest.mark.unit
    def test_protein_chain_case_insensitive(self):
        """Test ProteinChain accepts lowercase sequences."""
        chain = ProteinChain(["A"], "mvtp")
        assert chain.sequence == "mvtp"  # Should preserve original case

    @pytest.mark.unit
    def test_protein_chain_str(self, sample_protein_sequence):
        """Test ProteinChain string representation."""
        chain = ProteinChain(["A"], sample_protein_sequence, msa="test.a3m")
        str_repr = str(chain)
        assert "++ Protein chain:" in str_repr
        assert "MSA: 'test.a3m'" in str_repr

    @pytest.mark.unit
    def test_protein_chain_set_msa_path(self, sample_protein_sequence):
        """Test setting MSA path."""
        chain = ProteinChain(["A"], sample_protein_sequence)
        result = chain.set_msa_path("new_msa.a3m")
        assert chain.msa == "new_msa.a3m"
        assert result is chain  # Should return self

    @pytest.mark.unit
    def test_protein_chain_to_dict(self, sample_protein_sequence):
        """Test ProteinChain.to_dict()."""
        chain = ProteinChain(["A"], sample_protein_sequence, msa="test.a3m")
        result = chain.to_dict()

        expected = {
            "protein": {
                "id": "A",
                "sequence": sample_protein_sequence,
                "msa": "test.a3m",
            }
        }
        assert result == expected

    @pytest.mark.unit
    def test_protein_chain_to_dict_no_msa(self, sample_protein_sequence):
        """Test ProteinChain.to_dict() without MSA."""
        chain = ProteinChain(["A"], sample_protein_sequence)
        result = chain.to_dict()

        # MSA field should not be present if None
        assert "msa" not in result["protein"]

    @pytest.mark.unit
    def test_protein_chain_from_dict(self, sample_protein_sequence):
        """Test ProteinChain.from_dict()."""
        data = {"id": "A", "sequence": sample_protein_sequence, "msa": "test.a3m"}
        chain = ProteinChain.from_dict(data)

        assert chain.ids == ["A"]
        assert chain.sequence == sample_protein_sequence
        assert chain.msa == "test.a3m"


class TestDnaChain:
    """Test DnaChain class."""

    @pytest.mark.unit
    def test_dna_chain_creation(self, sample_dna_sequence):
        """Test basic DnaChain creation."""
        chain = DnaChain(["B"], sample_dna_sequence)
        assert chain.ids == ["B"]
        assert chain.sequence == sample_dna_sequence

    @pytest.mark.unit
    def test_dna_chain_invalid_sequence(self, invalid_dna_sequence):
        """Test DnaChain with invalid nucleotides."""
        with pytest.raises(
            ValueError, match="DNA sequence can only contain 'A', 'C', 'G', or 'T'"
        ):
            DnaChain(["B"], invalid_dna_sequence)

    @pytest.mark.unit
    def test_dna_chain_case_insensitive(self):
        """Test DnaChain accepts lowercase sequences."""
        chain = DnaChain(["B"], "atcg")
        assert chain.sequence == "atcg"  # Should preserve original case

    @pytest.mark.unit
    def test_dna_chain_str(self, sample_dna_sequence):
        """Test DnaChain string representation."""
        chain = DnaChain(["B"], sample_dna_sequence)
        str_repr = str(chain)
        assert "++ DNA chain:" in str_repr

    @pytest.mark.unit
    def test_dna_chain_to_dict(self, sample_dna_sequence):
        """Test DnaChain.to_dict()."""
        chain = DnaChain(["B"], sample_dna_sequence)
        result = chain.to_dict()

        expected = {"dna": {"id": "B", "sequence": sample_dna_sequence}}
        assert result == expected

    @pytest.mark.unit
    def test_dna_chain_from_dict(self, sample_dna_sequence):
        """Test DnaChain.from_dict()."""
        data = {"id": "B", "sequence": sample_dna_sequence}
        chain = DnaChain.from_dict(data)

        assert chain.ids == ["B"]
        assert chain.sequence == sample_dna_sequence


class TestRnaChain:
    """Test RnaChain class."""

    @pytest.mark.unit
    def test_rna_chain_creation(self, sample_rna_sequence):
        """Test basic RnaChain creation."""
        chain = RnaChain(["C"], sample_rna_sequence)
        assert chain.ids == ["C"]
        assert chain.sequence == sample_rna_sequence

    @pytest.mark.unit
    def test_rna_chain_invalid_sequence(self, invalid_rna_sequence):
        """Test RnaChain with invalid nucleotides."""
        with pytest.raises(
            ValueError, match="RNA sequence can only contain 'A', 'C', 'G', or 'U'"
        ):
            RnaChain(["C"], invalid_rna_sequence)

    @pytest.mark.unit
    def test_rna_chain_case_insensitive(self):
        """Test RnaChain accepts lowercase sequences."""
        chain = RnaChain(["C"], "aucg")
        assert chain.sequence == "aucg"  # Should preserve original case

    @pytest.mark.unit
    def test_rna_chain_str(self, sample_rna_sequence):
        """Test RnaChain string representation."""
        chain = RnaChain(["C"], sample_rna_sequence)
        str_repr = str(chain)
        assert "++ RNA chain:" in str_repr

    @pytest.mark.unit
    def test_rna_chain_to_dict(self, sample_rna_sequence):
        """Test RnaChain.to_dict()."""
        chain = RnaChain(["C"], sample_rna_sequence)
        result = chain.to_dict()

        expected = {"rna": {"id": "C", "sequence": sample_rna_sequence}}
        assert result == expected

    @pytest.mark.unit
    def test_rna_chain_from_dict(self, sample_rna_sequence):
        """Test RnaChain.from_dict()."""
        data = {"id": "C", "sequence": sample_rna_sequence}
        chain = RnaChain.from_dict(data)

        assert chain.ids == ["C"]
        assert chain.sequence == sample_rna_sequence


class TestLigand:
    """Test Ligand class."""

    @pytest.mark.unit
    def test_ligand_creation_smiles(self, sample_smiles):
        """Test Ligand creation with SMILES."""
        ligand = Ligand(["L"], smiles=sample_smiles)
        assert ligand.ids == ["L"]
        assert ligand.smiles == sample_smiles
        assert ligand.ccd is None

    @pytest.mark.unit
    def test_ligand_creation_ccd(self, sample_ccd):
        """Test Ligand creation with CCD."""
        ligand = Ligand(["L"], ccd=sample_ccd)
        assert ligand.ids == ["L"]
        assert ligand.ccd == sample_ccd
        assert ligand.smiles == ""

    @pytest.mark.unit
    def test_ligand_str_smiles(self, sample_smiles):
        """Test Ligand string representation with SMILES."""
        ligand = Ligand(["L"], smiles=sample_smiles)
        str_repr = str(ligand)
        assert "++ Ligand:" in str_repr
        assert f"SMILES: {sample_smiles}" in str_repr

    @pytest.mark.unit
    def test_ligand_str_ccd(self, sample_ccd):
        """Test Ligand string representation with CCD."""
        ligand = Ligand(["L"], ccd=sample_ccd)
        str_repr = str(ligand)
        assert "++ Ligand:" in str_repr
        assert f"CCD: {sample_ccd}" in str_repr

    @pytest.mark.unit
    def test_ligand_str_multiple_ids(self, sample_smiles):
        """Test Ligand string representation with multiple IDs."""
        ligand = Ligand(["L1", "L2"], smiles=sample_smiles)
        str_repr = str(ligand)
        assert "IDs: L1, L2" in str_repr

    @pytest.mark.unit
    def test_ligand_to_dict_smiles(self, sample_smiles):
        """Test Ligand.to_dict() with SMILES."""
        ligand = Ligand(["L"], smiles=sample_smiles)
        result = ligand.to_dict()

        expected = {"ligand": {"id": "L", "smiles": SingleQuoted(sample_smiles)}}
        assert result == expected
        assert isinstance(result["ligand"]["smiles"], SingleQuoted)

    @pytest.mark.unit
    def test_ligand_to_dict_ccd(self, sample_ccd):
        """Test Ligand.to_dict() with CCD."""
        ligand = Ligand(["L"], ccd=sample_ccd)
        result = ligand.to_dict()

        expected = {"ligand": {"id": "L", "ccd": sample_ccd}}
        assert result == expected

    @pytest.mark.unit
    def test_ligand_to_dict_multiple_ids(self, sample_smiles):
        """Test Ligand.to_dict() with multiple IDs."""
        ligand = Ligand(["L1", "L2"], smiles=sample_smiles)
        result = ligand.to_dict()

        assert isinstance(result["ligand"]["id"], FlowStyleList)
        assert result["ligand"]["id"] == ["L1", "L2"]

    @pytest.mark.unit
    def test_ligand_to_dict_no_smiles_no_ccd(self):
        """Test Ligand.to_dict() raises error with neither SMILES nor CCD."""
        ligand = Ligand(["L"])
        with pytest.raises(
            RuntimeError, match="Either a SMILES or a CCD code must be provided"
        ):
            ligand.to_dict()

    @pytest.mark.unit
    def test_ligand_from_dict_smiles(self, sample_smiles):
        """Test Ligand.from_dict() with SMILES."""
        data = {"id": "L", "smiles": sample_smiles}
        ligand = Ligand.from_dict(data)

        assert ligand.ids == ["L"]
        assert ligand.smiles == sample_smiles
        assert ligand.ccd is None

    @pytest.mark.unit
    def test_ligand_from_dict_ccd(self, sample_ccd):
        """Test Ligand.from_dict() with CCD."""
        data = {"id": "L", "ccd": sample_ccd}
        ligand = Ligand.from_dict(data)

        assert ligand.ids == ["L"]
        assert ligand.ccd == sample_ccd
        assert ligand.smiles == ""


class TestBond:
    """Test Bond constraint class."""

    @pytest.mark.unit
    def test_bond_creation(self, bond):
        """Test Bond creation."""
        assert bond.atom1 == ("A", 1, "CA")
        assert bond.atom2 == ("B", 2, "CB")

    @pytest.mark.unit
    def test_bond_str(self, bond):
        """Test Bond string representation."""
        str_repr = str(bond)
        assert "++ Bond between ('A', 1, 'CA') - ('B', 2, 'CB')" in str_repr

    @pytest.mark.unit
    def test_bond_to_dict(self, bond):
        """Test Bond.to_dict()."""
        result = bond.to_dict()

        expected = {
            "bond": {
                "atom1": FlowStyleList(("A", 1, "CA")),
                "atom2": FlowStyleList(("B", 2, "CB")),
            }
        }
        assert result == expected
        assert isinstance(result["bond"]["atom1"], FlowStyleList)
        assert isinstance(result["bond"]["atom2"], FlowStyleList)

    @pytest.mark.unit
    def test_bond_from_dict(self):
        """Test Bond.from_dict()."""
        data = {"atom1": ["A", 1, "CA"], "atom2": ["B", 2, "CB"]}
        bond = Bond.from_dict(data)

        assert bond.atom1 == ("A", 1, "CA")
        assert bond.atom2 == ("B", 2, "CB")


class TestContact:
    """Test Contact constraint class."""

    @pytest.mark.unit
    def test_contact_creation(self, contact):
        """Test Contact creation."""
        assert contact.token1 == ("A", 1)
        assert contact.token2 == ("B", 2)
        assert contact.max_distance == 6.0

    @pytest.mark.unit
    def test_contact_creation_default_distance(self):
        """Test contact creation with default distance."""
        contact = Contact(("A", 1), ("B", 2))
        assert contact.max_distance == 6.0

    @pytest.mark.unit
    def test_contact_str(self, contact):
        """Test Contact string representation."""
        str_repr = str(contact)
        assert (
            "++ Contact between ('A', 1) - ('B', 2), max distance: 6.00 Å" in str_repr
        )

    @pytest.mark.unit
    def test_contact_to_dict(self, contact):
        """Test Contact.to_dict()."""
        result = contact.to_dict()

        expected = {
            "contact": {
                "token1": FlowStyleList(("A", 1)),
                "token2": FlowStyleList(("B", 2)),
                "max_distance": 6.0,
            }
        }
        assert result == expected

    @pytest.mark.unit
    def test_contact_from_dict(self):
        """Test Contact.from_dict()."""
        data = {"token1": ["A", 1], "token2": ["B", 2], "max_distance": 6.0}
        contact = Contact.from_dict(data)

        assert contact.token1 == ("A", 1)
        assert contact.token2 == ("B", 2)
        assert contact.max_distance == 6.0


class TestPocket:
    """Test Pocket constraint class."""

    @pytest.mark.unit
    def test_pocket_creation(self, pocket):
        """Test Pocket creation."""
        assert pocket.binder == "L"
        assert pocket.contacts == [("A", 5), ("A", 10)]
        assert pocket.max_distance == 5.0

    @pytest.mark.unit
    def test_pocket_creation_defaults(self):
        """Test pocket creation with defaults."""
        pocket = Pocket("L")
        assert pocket.binder == "L"
        assert pocket.contacts == []
        assert pocket.max_distance == 6.0

    @pytest.mark.unit
    def test_pocket_str(self, pocket):
        """Test Pocket string representation."""
        str_repr = str(pocket)
        assert "++ Pocket for binder: L" in str_repr
        assert "Max distance: 5.00 Å" in str_repr
        assert "Contacts:" in str_repr
        assert "- A 5" in str_repr
        assert "- A 10" in str_repr

    @pytest.mark.unit
    def test_pocket_str_no_contacts(self):
        """Test Pocket string representation with no contacts."""
        pocket = Pocket("L")
        str_repr = str(pocket)
        assert "No contacts defined." in str_repr

    @pytest.mark.unit
    def test_pocket_add_contact_token(self, pocket):
        """Test adding contact tokens."""
        initial_contacts = len(pocket.contacts)
        pocket.add_contact_token("B", 15)

        assert len(pocket.contacts) == initial_contacts + 1
        assert ("B", 15) in pocket.contacts

    @pytest.mark.unit
    def test_pocket_to_dict(self, pocket):
        """Test Pocket.to_dict()."""
        result = pocket.to_dict()

        expected = {
            "pocket": {
                "binder": "L",
                "contacts": FlowStyleList([("A", 5), ("A", 10)]),
                "max_distance": 5.0,
            }
        }
        assert result == expected

    @pytest.mark.unit
    def test_pocket_from_dict(self):
        """Test Pocket.from_dict()."""
        data = {"binder": "L", "contacts": [["A", 5], ["A", 10]], "max_distance": 5.0}
        pocket = Pocket.from_dict(data)

        assert pocket.binder == "L"
        assert pocket.contacts == [["A", 5], ["A", 10]]
        assert pocket.max_distance == 5.0


class TestTemplate:
    """Test Template class."""

    @pytest.mark.unit
    def test_template_creation(self, template):
        """Test Template creation."""
        assert template.cif == "test.cif"
        assert template.chain_id == ["A"]
        assert template.template_id == ["X"]

    @pytest.mark.unit
    def test_template_creation_defaults(self):
        """Test Template creation with defaults."""
        template = Template("test.cif")
        assert template.chain_id == []
        assert template.template_id == []

    @pytest.mark.unit
    def test_template_str(self, template):
        """Test Template string representation."""
        str_repr = str(template)
        assert "'test.cif'" in str_repr
        assert "Template ID(s): X" in str_repr
        assert "For chain ID(s): A" in str_repr

    @pytest.mark.unit
    def test_template_to_dict_single_ids(self, template):
        """Test Template.to_dict() with single IDs."""
        result = template.to_dict()

        expected = {
            "cif": "test.cif",
            "chain_id": "A",  # Single ID as string
            "template_id": "X",  # Single ID as string
        }
        assert result == expected

    @pytest.mark.unit
    def test_template_to_dict_multiple_ids(self):
        """Test Template.to_dict() with multiple IDs."""
        template = Template("test.cif", chain_id=["A", "B"], template_id=["X", "Y"])
        result = template.to_dict()

        expected = {
            "cif": "test.cif",
            "chain_id": FlowStyleList(["A", "B"]),
            "template_id": FlowStyleList(["X", "Y"]),
        }
        assert result == expected
        assert isinstance(result["chain_id"], FlowStyleList)
        assert isinstance(result["template_id"], FlowStyleList)

    @pytest.mark.unit
    def test_template_to_dict_no_optional_fields(self):
        """Test Template.to_dict() without optional fields."""
        template = Template("test.cif")
        result = template.to_dict()

        assert result == {"cif": "test.cif"}
        assert "chain_id" not in result
        assert "template_id" not in result

    @pytest.mark.unit
    def test_template_from_dict(self):
        """Test Template.from_dict()."""
        data = {"cif": "test.cif", "chain_id": "A", "template_id": "X"}
        template = Template.from_dict(data)

        assert template.cif == "test.cif"
        assert template.chain_id == ["A"]
        assert template.template_id == ["X"]

    @pytest.mark.unit
    def test_template_from_dict_multiple_ids(self):
        """Test Template.from_dict() with multiple IDs."""
        data = {"cif": "test.cif", "chain_id": ["A", "B"], "template_id": ["X", "Y"]}
        template = Template.from_dict(data)

        assert template.chain_id == ["A", "B"]
        assert template.template_id == ["X", "Y"]


class TestAffinity:
    """Test Affinity property class."""

    @pytest.mark.unit
    def test_affinity_creation(self, affinity):
        """Test Affinity creation."""
        assert affinity.binder == "L"

    @pytest.mark.unit
    def test_affinity_str(self, affinity):
        """Test Affinity string representation."""
        str_repr = str(affinity)
        assert "++ Affinity for ligand: L" in str_repr

    @pytest.mark.unit
    def test_affinity_to_dict(self, affinity):
        """Test Affinity.to_dict()."""
        result = affinity.to_dict()

        expected = {"affinity": {"binder": "L"}}
        assert result == expected

    @pytest.mark.unit
    def test_affinity_from_dict(self):
        """Test Affinity.from_dict()."""
        data = {"binder": "L"}
        affinity = Affinity.from_dict(data)

        assert affinity.binder == "L"

    @pytest.mark.unit
    def test_affinity_from_dict_empty(self):
        """Test Affinity.from_dict() with empty data."""
        affinity = Affinity.from_dict({})
        assert affinity.binder == ""

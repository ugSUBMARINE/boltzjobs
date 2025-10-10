"""Tests for boltzjobs.jobs module."""

import pytest
import warnings
from pathlib import Path

from boltzjobs import Job
from boltzjobs.components import (
    ProteinChain, DnaChain, RnaChain, Ligand,
    Bond, Contact, Pocket, Template, Affinity
)


class TestJobCreation:
    """Test Job class creation and basic properties."""
    
    @pytest.mark.unit
    def test_job_default_creation(self):
        """Test Job creation with defaults."""
        job = Job()
        assert job.name == "boltz_job"
        assert job.version == 1
        assert job.sequences == []
        assert job.constraints == []
        assert job.properties == []
        assert job.templates == []
    
    @pytest.mark.unit
    def test_job_custom_creation(self):
        """Test Job creation with custom parameters."""
        job = Job("custom_job", version=2)
        assert job.name == "custom_job"
        assert job.version == 2
    
    @pytest.mark.unit
    def test_job_post_init(self):
        """Test Job post-initialization creates chain ID generator."""
        job = Job()
        # Should have initialized chain ID generator
        assert hasattr(job, '_chain_ids')
        # Should be able to get chain IDs
        ids = job._get_ids(None, 3)
        assert ids == ["A", "B", "C"]


class TestJobStringRepresentation:
    """Test Job string representation."""
    
    @pytest.mark.unit
    def test_job_str_empty(self):
        """Test Job string representation when empty."""
        job = Job("test_job", version=2)
        str_repr = str(job)
        assert "Job name: test_job" in str_repr
        assert "Version: 2" in str_repr
        assert "No sequences." in str_repr
    
    @pytest.mark.unit
    def test_job_str_with_sequences(self, basic_job):
        """Test Job string representation with sequences."""
        str_repr = str(basic_job)
        assert "Job name: test_job" in str_repr
        assert "2 Sequence(s):" in str_repr
        assert "++ Protein chain:" in str_repr
        assert "++ Ligand:" in str_repr
    
    @pytest.mark.unit
    def test_job_str_complex(self, complex_job):
        """Test Job string representation with all components."""
        str_repr = str(complex_job)
        assert "Job name: complex_test_job" in str_repr
        assert "3 Sequence(s):" in str_repr
        assert "Constraints:" in str_repr
        assert "Templates:" in str_repr
        assert "Properties:" in str_repr


class TestJobIdManagement:
    """Test Job ID management functionality."""
    
    @pytest.mark.unit
    def test_get_current_ids(self, basic_job):
        """Test getting current chain IDs."""
        ids = basic_job._get_current_ids()
        assert "A" in ids
        assert "L" in ids
        assert len(ids) == 2
    
    @pytest.mark.unit
    def test_check_ids_no_duplicates(self, basic_job):
        """Test ID validation with no duplicates."""
        # Should not raise
        basic_job._check_ids()
    
    @pytest.mark.unit
    def test_check_ids_with_duplicates(self, sample_protein_sequence):
        """Test ID validation with duplicates raises error."""
        job = Job()
        # Add two sequences with same ID manually
        job.sequences = [
            ProteinChain(["A"], sample_protein_sequence),
            ProteinChain(["A"], sample_protein_sequence)
        ]
        
        with pytest.raises(ValueError, match="Duplicate chain IDs found"):
            job._check_ids()
    
    @pytest.mark.unit
    def test_get_ids_none(self):
        """Test getting IDs when none provided."""
        job = Job()
        ids = job._get_ids(None, 3)
        assert ids == ["A", "B", "C"]
    
    @pytest.mark.unit
    def test_get_ids_string(self):
        """Test getting IDs from string."""
        job = Job()
        ids = job._get_ids("X", 1)
        assert ids == ["X"]
    
    @pytest.mark.unit
    def test_get_ids_list(self):
        """Test getting IDs from list."""
        job = Job()
        ids = job._get_ids(["X", "Y"], 2)
        assert ids == ["X", "Y"]
    
    @pytest.mark.unit
    def test_get_ids_invalid_type(self):
        """Test getting IDs with invalid type."""
        job = Job()
        with pytest.raises(TypeError, match="IDs must be a string or a list of strings"):
            job._get_ids(123, 1)
    
    @pytest.mark.unit
    def test_get_ids_empty(self):
        """Test getting IDs with empty list."""
        job = Job()
        with pytest.raises(ValueError, match="Number of chains or ligands must be greater than zero"):
            job._get_ids([], 0)


class TestJobAddProteinChain:
    """Test Job.add_protein_chain() method."""
    
    @pytest.mark.unit
    def test_add_protein_chain_basic(self, sample_protein_sequence):
        """Test basic protein chain addition."""
        job = Job()
        chain = job.add_protein_chain(sample_protein_sequence)
        
        assert len(job.sequences) == 1
        assert isinstance(job.sequences[0], ProteinChain)
        assert job.sequences[0].sequence == sample_protein_sequence
        assert job.sequences[0].ids == ["A"]
        assert chain is job.sequences[0]
    
    @pytest.mark.unit
    def test_add_protein_chain_with_custom_id(self, sample_protein_sequence):
        """Test protein chain addition with custom ID."""
        job = Job()
        chain = job.add_protein_chain(sample_protein_sequence, ids="P1")
        
        assert chain.ids == ["P1"]
    
    @pytest.mark.unit
    def test_add_protein_chain_multiple_copies(self, sample_protein_sequence):
        """Test protein chain addition with multiple copies."""
        job = Job()
        chain = job.add_protein_chain(sample_protein_sequence, count=3)
        
        assert chain.ids == ["A", "B", "C"]
    
    @pytest.mark.unit
    def test_add_protein_chain_with_msa(self, sample_protein_sequence):
        """Test protein chain addition with MSA."""
        job = Job()
        chain = job.add_protein_chain(sample_protein_sequence, msa="test.a3m")
        
        assert chain.msa == "test.a3m"
    
    @pytest.mark.unit
    def test_add_protein_chain_cyclic(self, sample_protein_sequence):
        """Test protein chain addition with cyclic flag."""
        job = Job()
        chain = job.add_protein_chain(sample_protein_sequence, cyclic=True)
        
        assert chain.cyclic is True


class TestJobAddDnaChain:
    """Test Job.add_dna_chain() method."""
    
    @pytest.mark.unit
    def test_add_dna_chain_basic(self, sample_dna_sequence):
        """Test basic DNA chain addition."""
        job = Job()
        chain = job.add_dna_chain(sample_dna_sequence)
        
        assert len(job.sequences) == 1
        assert isinstance(job.sequences[0], DnaChain)
        assert job.sequences[0].sequence == sample_dna_sequence
        assert job.sequences[0].ids == ["A"]
    
    @pytest.mark.unit
    def test_add_dna_chain_with_custom_id(self, sample_dna_sequence):
        """Test DNA chain addition with custom ID."""
        job = Job()
        chain = job.add_dna_chain(sample_dna_sequence, ids="D1")
        
        assert chain.ids == ["D1"]


class TestJobAddRnaChain:
    """Test Job.add_rna_chain() method."""
    
    @pytest.mark.unit
    def test_add_rna_chain_basic(self, sample_rna_sequence):
        """Test basic RNA chain addition."""
        job = Job()
        chain = job.add_rna_chain(sample_rna_sequence)
        
        assert len(job.sequences) == 1
        assert isinstance(job.sequences[0], RnaChain)
        assert job.sequences[0].sequence == sample_rna_sequence
        assert job.sequences[0].ids == ["A"]
    
    @pytest.mark.unit
    def test_add_rna_chain_with_custom_id(self, sample_rna_sequence):
        """Test RNA chain addition with custom ID."""
        job = Job()
        chain = job.add_rna_chain(sample_rna_sequence, ids="R1")
        
        assert chain.ids == ["R1"]


class TestJobAddLigand:
    """Test Job.add_ligand() method."""
    
    @pytest.mark.unit
    def test_add_ligand_smiles(self, sample_smiles):
        """Test ligand addition with SMILES."""
        job = Job()
        ligand = job.add_ligand(smiles=sample_smiles)
        
        assert len(job.sequences) == 1
        assert isinstance(job.sequences[0], Ligand)
        assert job.sequences[0].smiles == sample_smiles
        assert job.sequences[0].ids == ["A"]
    
    @pytest.mark.unit
    def test_add_ligand_ccd(self, sample_ccd):
        """Test ligand addition with CCD."""
        job = Job()
        ligand = job.add_ligand(ccd=sample_ccd)
        
        assert job.sequences[0].ccd == sample_ccd
    
    @pytest.mark.unit
    def test_add_ligand_both_smiles_and_ccd(self, sample_smiles, sample_ccd):
        """Test ligand addition with both SMILES and CCD warns and uses SMILES."""
        job = Job()
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            ligand = job.add_ligand(smiles=sample_smiles, ccd=sample_ccd)
            
            assert len(w) == 1
            assert "mutually exclusive" in str(w[0].message)
            assert ligand.smiles == sample_smiles
            assert ligand.ccd is None
    
    @pytest.mark.unit
    def test_add_ligand_neither_smiles_nor_ccd(self):
        """Test ligand addition without SMILES or CCD raises error."""
        job = Job()
        with pytest.raises(ValueError, match="Either CCD codes or SMILES string must be provided"):
            job.add_ligand()
    
    @pytest.mark.unit
    def test_add_ligand_multiple_copies(self, sample_smiles):
        """Test ligand addition with multiple copies."""
        job = Job()
        ligand = job.add_ligand(smiles=sample_smiles, count=2)
        
        assert ligand.ids == ["A", "B"]


class TestJobAddConstraints:
    """Test Job constraint addition methods."""
    
    @pytest.mark.unit
    def test_add_bond(self, basic_job):
        """Test adding bond constraint."""
        basic_job.add_bond("A", 1, "CA", "L", 1, "O1")
        
        assert len(basic_job.constraints) == 1
        assert isinstance(basic_job.constraints[0], Bond)
        assert basic_job.constraints[0].atom1 == ("A", 1, "CA")
        assert basic_job.constraints[0].atom2 == ("L", 1, "O1")
    
    @pytest.mark.unit
    def test_add_bond_invalid_chain_ids(self):
        """Test adding bond with invalid chain IDs."""
        job = Job()
        with pytest.raises(ValueError, match="Both chain IDs must be defined"):
            job.add_bond("X", 1, "CA", "Y", 1, "O1")
    
    @pytest.mark.unit
    def test_add_disulfide_bond(self, sample_protein_sequence):
        """Test adding disulfide bond constraint."""
        # Create sequence with cysteines
        cys_seq = sample_protein_sequence.replace("M", "C").replace("V", "C")
        job = Job()
        job.add_protein_chain(cys_seq, ids="A")
        job.add_protein_chain(cys_seq, ids="B")
        
        job.add_disulfide_bond("A", 1, "B", 2)
        
        assert len(job.constraints) == 1
        assert isinstance(job.constraints[0], Bond)
        assert job.constraints[0].atom1[2] == "SG"
        assert job.constraints[0].atom2[2] == "SG"
    
    @pytest.mark.unit
    def test_add_disulfide_bond_non_cysteine(self, sample_protein_sequence):
        """Test adding disulfide bond to non-cysteine raises error."""
        job = Job()
        job.add_protein_chain(sample_protein_sequence, ids="A")
        job.add_protein_chain(sample_protein_sequence, ids="B")
        
        with pytest.raises(ValueError, match="is not a cysteine"):
            job.add_disulfide_bond("A", 1, "B", 1)  # M and V are not cysteines
    
    @pytest.mark.unit
    def test_add_contact(self, basic_job):
        """Test adding contact constraint."""
        basic_job.add_contact("A", 1, "L", 1, max_distance=6.0)
        
        assert len(basic_job.constraints) == 1
        assert isinstance(basic_job.constraints[0], Contact)
        assert basic_job.constraints[0].max_distance == 6.0
    
    @pytest.mark.unit
    def test_add_contact_invalid_chain_ids(self):
        """Test adding contact with invalid chain IDs."""
        job = Job()
        with pytest.raises(ValueError, match="Both chain IDs must be defined"):
            job.add_contact("X", 1, "Y", 1)
    
    @pytest.mark.unit
    def test_add_pocket(self, basic_job):
        """Test adding pocket constraint."""
        pocket = basic_job.add_pocket("L", contact_tokens=[("A", 5)], max_distance=6.0)
        
        assert len(basic_job.constraints) == 1
        assert isinstance(basic_job.constraints[0], Pocket)
        assert basic_job.constraints[0].binder == "L"
        assert basic_job.constraints[0].max_distance == 6.0
        assert pocket is basic_job.constraints[0]
    
    @pytest.mark.unit
    def test_add_pocket_allows_all_sequence_types(self, sample_protein_sequence):
        """Test adding pocket now works for all sequence types."""
        job = Job()
        job.add_protein_chain(sample_protein_sequence, ids="A")
        
        # Should now work for protein chains (changed behavior)
        pocket = job.add_pocket("A")
        assert pocket.binder == "A"
    
    @pytest.mark.unit
    def test_add_pocket_non_existent_binder(self, basic_job):
        """Test adding pocket with non-existent binder raises error."""
        with pytest.raises(ValueError, match="No sequence with id 'X' found"):
            basic_job.add_pocket("X")


class TestJobAddTemplate:
    """Test Job.add_template() method."""
    
    @pytest.mark.unit
    def test_add_template_basic(self):
        """Test basic template addition."""
        job = Job()
        job.add_template("test.cif")
        
        assert len(job.templates) == 1
        assert isinstance(job.templates[0], Template)
        assert job.templates[0].cif == "test.cif"
        assert job.templates[0].chain_id == []
        assert job.templates[0].template_id == []
    
    @pytest.mark.unit
    def test_add_template_with_chain_id(self):
        """Test template addition with chain ID."""
        job = Job()
        job.add_template("test.cif", chain_id="A")
        
        assert job.templates[0].chain_id == ["A"]
    
    @pytest.mark.unit
    def test_add_template_with_multiple_chain_ids(self):
        """Test template addition with multiple chain IDs."""
        job = Job()
        job.add_template("test.cif", chain_id=["A", "B"], template_id=["X", "Y"])
        
        assert job.templates[0].chain_id == ["A", "B"]
        assert job.templates[0].template_id == ["X", "Y"]
    
    @pytest.mark.unit
    def test_add_template_empty_cif(self):
        """Test template addition with empty CIF raises error."""
        job = Job()
        with pytest.raises(ValueError, match="Either 'cif' or 'pdb' file path must be provided"):
            job.add_template("")


class TestJobRequestAffinity:
    """Test Job.request_affinity() method."""
    
    @pytest.mark.unit
    def test_request_affinity(self, basic_job):
        """Test requesting affinity for ligand."""
        basic_job.request_affinity("L")
        
        assert len(basic_job.properties) == 1
        assert isinstance(basic_job.properties[0], Affinity)
        assert basic_job.properties[0].binder == "L"
    
    @pytest.mark.unit
    def test_request_affinity_non_ligand(self, sample_protein_sequence):
        """Test requesting affinity for non-ligand raises error."""
        job = Job()
        job.add_protein_chain(sample_protein_sequence, ids="A")
        
        with pytest.raises(ValueError, match="Affinity can only be estimated for ligands"):
            job.request_affinity("A")
    
    @pytest.mark.unit
    def test_request_affinity_non_existent_ligand(self, basic_job):
        """Test requesting affinity for non-existent ligand raises error."""
        with pytest.raises(ValueError, match="No ligand with id 'X' found"):
            basic_job.request_affinity("X")


class TestJobToDictYaml:
    """Test Job.to_dict() and YAML functionality."""
    
    @pytest.mark.unit
    def test_to_dict_basic(self, basic_job):
        """Test basic to_dict functionality."""
        result = basic_job.to_dict()
        
        assert "version" in result
        assert result["version"] == 1
        assert "sequences" in result
        assert len(result["sequences"]) == 2
        
        # Check protein chain
        protein_seq = None
        ligand_seq = None
        for seq in result["sequences"]:
            if "protein" in seq:
                protein_seq = seq["protein"]
            elif "ligand" in seq:
                ligand_seq = seq["ligand"]
        
        assert protein_seq is not None
        assert ligand_seq is not None
        assert protein_seq["id"] == "A"
        assert ligand_seq["id"] == "L"
    
    @pytest.mark.unit
    def test_to_dict_complex(self, complex_job):
        """Test to_dict with all components."""
        result = complex_job.to_dict()
        
        assert "version" in result
        assert "sequences" in result
        assert "constraints" in result
        assert "templates" in result
        assert "properties" in result
        
        assert len(result["sequences"]) == 3
        assert len(result["constraints"]) == 3
        assert len(result["templates"]) == 1
        assert len(result["properties"]) == 1
    
    @pytest.mark.unit
    def test_to_dict_empty_sequences_error(self):
        """Test to_dict with empty sequences raises error."""
        job = Job()
        with pytest.raises(ValueError, match="Empty list of sequences"):
            job.to_dict()
    
    @pytest.mark.unit
    def test_to_dict_duplicate_ids_error(self, sample_protein_sequence):
        """Test to_dict with duplicate IDs raises error."""
        job = Job()
        job.sequences = [
            ProteinChain(["A"], sample_protein_sequence),
            ProteinChain(["A"], sample_protein_sequence)
        ]
        
        with pytest.raises(ValueError, match="Duplicate chain IDs found"):
            job.to_dict()
    
    @pytest.mark.unit
    def test_write_yaml(self, basic_job, temp_yaml_file):
        """Test writing job to YAML file."""
        basic_job.write_yaml(temp_yaml_file)
        
        # File should exist and be readable
        assert Path(temp_yaml_file).exists()
        
        # Should be able to read back
        with open(temp_yaml_file) as f:
            content = f.read()
            assert "version: 1" in content
            assert "sequences:" in content


class TestJobFromYaml:
    """Test Job.from_yaml() functionality."""
    
    @pytest.mark.unit
    def test_from_yaml_basic(self, temp_yaml_file, sample_yaml_content):
        """Test loading job from YAML."""
        # Write sample YAML content
        with open(temp_yaml_file, 'w') as f:
            f.write(sample_yaml_content)
        
        job = Job.from_yaml(temp_yaml_file, "loaded_job")
        
        assert job.name == "loaded_job"
        assert job.version == 1
        assert len(job.sequences) == 2
        assert len(job.constraints) == 1
        assert len(job.properties) == 1
        
        # Check sequences
        protein_chain = None
        ligand = None
        for seq in job.sequences:
            if isinstance(seq, ProteinChain):
                protein_chain = seq
            elif isinstance(seq, Ligand):
                ligand = seq
        
        assert protein_chain is not None
        assert ligand is not None
        assert protein_chain.ids == ["A"]
        assert ligand.ids == ["L"]
    
    @pytest.mark.unit
    def test_from_yaml_version_default(self, temp_yaml_file):
        """Test loading job from YAML with default version."""
        yaml_content = """sequences:
  - protein:
      id: A
      sequence: MVTP
"""
        with open(temp_yaml_file, 'w') as f:
            f.write(yaml_content)
        
        job = Job.from_yaml(temp_yaml_file)
        assert job.version == 1  # Should default to 1
    
    @pytest.mark.unit
    def test_from_yaml_empty_sections(self, temp_yaml_file):
        """Test loading job from YAML with only sequences."""
        yaml_content = """version: 1
sequences:
  - protein:
      id: A
      sequence: MVTP
"""
        with open(temp_yaml_file, 'w') as f:
            f.write(yaml_content)
        
        job = Job.from_yaml(temp_yaml_file)
        assert len(job.sequences) == 1
        assert len(job.constraints) == 0
        assert len(job.templates) == 0
        assert len(job.properties) == 0


class TestJobRoundTrip:
    """Test YAML round-trip functionality."""
    
    @pytest.mark.unit
    def test_yaml_round_trip_basic(self, basic_job, temp_yaml_file):
        """Test basic YAML round-trip."""
        # Write to YAML
        basic_job.write_yaml(temp_yaml_file)
        
        # Read back
        loaded_job = Job.from_yaml(temp_yaml_file, basic_job.name)
        
        # Should have same basic structure
        assert loaded_job.name == basic_job.name
        assert loaded_job.version == basic_job.version
        assert len(loaded_job.sequences) == len(basic_job.sequences)
    
    @pytest.mark.unit  
    def test_yaml_round_trip_complex(self, complex_job, temp_yaml_file):
        """Test complex YAML round-trip."""
        # Write to YAML
        complex_job.write_yaml(temp_yaml_file)
        
        # Read back
        loaded_job = Job.from_yaml(temp_yaml_file, complex_job.name)
        
        # Should have same structure
        assert loaded_job.name == complex_job.name
        assert len(loaded_job.sequences) == len(complex_job.sequences)
        assert len(loaded_job.constraints) == len(complex_job.constraints)
        assert len(loaded_job.templates) == len(complex_job.templates)
        assert len(loaded_job.properties) == len(complex_job.properties)
        
        # Check specific components
        assert isinstance(loaded_job.sequences[0], ProteinChain)
        assert isinstance(loaded_job.sequences[1], DnaChain)
        assert isinstance(loaded_job.sequences[2], Ligand)
        
        # Check constraint types
        constraint_types = [type(c).__name__ for c in loaded_job.constraints]
        assert "Bond" in constraint_types
        assert "Contact" in constraint_types
        assert "Pocket" in constraint_types
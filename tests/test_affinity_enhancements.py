"""Tests for Enhanced Affinity Computation (Section 4 of TODO.md)."""

import pytest
import warnings

from boltzjobs import Job
from boltzjobs.components import Affinity


class TestAffinityValidation:
    """Test enhanced affinity validation functionality."""

    @pytest.fixture
    def basic_protein_ligand_job(self):
        """Create a basic job with protein and ligand for testing."""
        job = Job("affinity_test")
        job.add_protein_chain("MVTPKVLH", ids=["A"])
        job.add_ligand(smiles="CCO", ids=["LIG"])
        return job

    @pytest.fixture
    def dna_rna_ligand_job(self):
        """Create a job with DNA, RNA, and ligand (no protein)."""
        job = Job("non_protein_test")
        job.add_dna_chain("ATCG", ids=["DNA"])
        job.add_rna_chain("AUCG", ids=["RNA"]) 
        job.add_ligand(smiles="CCO", ids=["LIG"])
        return job

    @pytest.fixture
    def ligand_only_job(self):
        """Create a job with only ligands (no protein targets)."""
        job = Job("ligand_only_test")
        job.add_ligand(smiles="CCO", ids=["LIG1"])
        job.add_ligand(ccd="ATP", ids=["LIG2"])
        return job

    @pytest.fixture
    def mixed_targets_job(self):
        """Create a job with protein, DNA, RNA, and ligand."""
        job = Job("mixed_targets_test")
        job.add_protein_chain("MVTP", ids=["PROT"])
        job.add_dna_chain("ATCG", ids=["DNA"])
        job.add_rna_chain("AUCG", ids=["RNA"])
        job.add_ligand(smiles="CCO", ids=["LIG"])
        return job

    @pytest.mark.unit
    def test_request_affinity_basic(self, basic_protein_ligand_job):
        """Test basic affinity request functionality."""
        basic_protein_ligand_job.request_affinity("LIG")
        
        assert len(basic_protein_ligand_job.properties) == 1
        assert isinstance(basic_protein_ligand_job.properties[0], Affinity)
        assert basic_protein_ligand_job.properties[0].binder == "LIG"

    @pytest.mark.unit
    def test_request_affinity_non_existent_ligand(self, basic_protein_ligand_job):
        """Test affinity request for non-existent ligand raises error."""
        with pytest.raises(ValueError, match="No ligand with id 'NONEXISTENT' found"):
            basic_protein_ligand_job.request_affinity("NONEXISTENT")

    @pytest.mark.unit
    def test_request_affinity_for_protein(self, basic_protein_ligand_job):
        """Test affinity request for protein chain raises error."""
        with pytest.raises(ValueError, match="No ligand with id 'A' found"):
            basic_protein_ligand_job.request_affinity("A")

    @pytest.mark.unit
    def test_request_affinity_duplicate_error(self, basic_protein_ligand_job):
        """Test that requesting affinity twice raises error."""
        # First request should succeed
        basic_protein_ligand_job.request_affinity("LIG")
        assert len(basic_protein_ligand_job.properties) == 1
        
        # Second request should fail
        with pytest.raises(ValueError, match="Only one affinity computation allowed per job"):
            basic_protein_ligand_job.request_affinity("LIG")

    @pytest.mark.unit
    def test_request_affinity_duplicate_different_ligand(self):
        """Test that requesting affinity for different ligands fails."""
        job = Job("multi_ligand_test")
        job.add_protein_chain("MVTP", ids=["A"])
        job.add_ligand(smiles="CCO", ids=["LIG1"])
        job.add_ligand(ccd="ATP", ids=["LIG2"])
        
        # First request should succeed
        job.request_affinity("LIG1")
        assert len(job.properties) == 1
        assert job.properties[0].binder == "LIG1"
        
        # Second request for different ligand should fail
        with pytest.raises(ValueError, match="Affinity already requested for: LIG1"):
            job.request_affinity("LIG2")

    @pytest.mark.unit
    def test_request_affinity_with_protein_targets_no_warning(self, basic_protein_ligand_job):
        """Test affinity request with protein targets doesn't issue warnings."""
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            basic_protein_ligand_job.request_affinity("LIG")
            
            # Should not have any warnings for protein targets
            assert len(w) == 0

    @pytest.mark.unit
    def test_request_affinity_no_protein_targets_warning(self, ligand_only_job):
        """Test affinity request without protein targets issues warning."""
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            ligand_only_job.request_affinity("LIG1")
            
            # Should have warning about no protein targets
            assert len(w) == 1
            assert "No protein chains found" in str(w[0].message)
            assert "may produce unreliable results" in str(w[0].message)
            assert w[0].category == UserWarning

    @pytest.mark.unit
    def test_request_affinity_dna_rna_targets_warning(self, dna_rna_ligand_job):
        """Test affinity request with DNA/RNA targets issues warnings."""
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            dna_rna_ligand_job.request_affinity("LIG")
            
            # Should have two warnings: no protein + non-protein targets
            assert len(w) == 2
            
            # Check no protein warning
            no_protein_warning = next((warning for warning in w 
                                     if "No protein chains found" in str(warning.message)), None)
            assert no_protein_warning is not None
            
            # Check non-protein targets warning
            non_protein_warning = next((warning for warning in w 
                                      if "Non-protein targets detected" in str(warning.message)), None)
            assert non_protein_warning is not None
            assert "DnaChain, RnaChain" in str(non_protein_warning.message)
            assert "may produce unreliable results" in str(non_protein_warning.message)

    @pytest.mark.unit
    def test_request_affinity_mixed_targets_warning(self, mixed_targets_job):
        """Test affinity request with mixed targets issues appropriate warning."""
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            mixed_targets_job.request_affinity("LIG")
            
            # Should have one warning about non-protein targets (but not about missing protein)
            assert len(w) == 1
            warning = w[0]
            assert "Non-protein targets detected" in str(warning.message)
            assert "DnaChain, RnaChain" in str(warning.message)
            assert w[0].category == UserWarning

    @pytest.mark.unit
    def test_request_affinity_preserves_existing_properties(self):
        """Test that requesting affinity preserves other existing properties."""
        job = Job("preserve_test")
        job.add_protein_chain("MVTP", ids=["A"])
        job.add_ligand(smiles="CCO", ids=["LIG"])
        
        # Manually add a non-affinity property (simulating future property types)
        # Since we only have Affinity now, we'll test the general case
        job.request_affinity("LIG")
        
        # Verify the affinity was added correctly
        assert len(job.properties) == 1
        assert isinstance(job.properties[0], Affinity)
        assert job.properties[0].binder == "LIG"

    @pytest.mark.unit
    def test_request_affinity_multiple_ligands_same_id(self):
        """Test affinity request works with ligands that have multiple IDs."""
        job = Job("multi_id_test")
        job.add_protein_chain("MVTP", ids=["A"])
        job.add_ligand(smiles="CCO", ids=["LIG", "ETHANOL"])
        
        # Should work with either ID
        job.request_affinity("ETHANOL")
        assert len(job.properties) == 1
        assert job.properties[0].binder == "ETHANOL"

    @pytest.mark.unit
    def test_request_affinity_case_sensitivity(self, basic_protein_ligand_job):
        """Test that affinity request is case sensitive."""
        with pytest.raises(ValueError, match="No ligand with id 'lig' found"):
            basic_protein_ligand_job.request_affinity("lig")  # lowercase

    @pytest.mark.unit
    def test_request_affinity_empty_job(self):
        """Test affinity request on empty job raises appropriate error."""
        job = Job("empty_test")
        
        with pytest.raises(ValueError, match="No ligand with id 'LIG' found"):
            job.request_affinity("LIG")


class TestAffinityIntegration:
    """Test affinity integration with job workflows."""

    @pytest.mark.unit
    def test_affinity_yaml_serialization(self):
        """Test that affinity request serializes correctly to YAML."""
        job = Job("yaml_test")
        job.add_protein_chain("MVTP", ids=["A"])
        job.add_ligand(smiles="CCO", ids=["LIG"])
        job.request_affinity("LIG")
        
        yaml_dict = job.to_dict()
        
        assert "properties" in yaml_dict
        assert len(yaml_dict["properties"]) == 1
        assert "affinity" in yaml_dict["properties"][0]
        assert yaml_dict["properties"][0]["affinity"]["binder"] == "LIG"

    @pytest.mark.unit
    def test_affinity_with_constraints(self):
        """Test affinity works alongside constraints."""
        job = Job("constraints_test")
        job.add_protein_chain("MVTPKVLHC", ids=["A"])
        job.add_ligand(smiles="CCO", ids=["LIG"])
        
        # Add constraints
        job.add_pocket("LIG")
        job.add_contact("A", 1, "LIG", 1)
        
        # Request affinity
        job.request_affinity("LIG")
        
        # Verify all components are present
        assert len(job.constraints) == 2
        assert len(job.properties) == 1
        assert job.properties[0].binder == "LIG"

    @pytest.mark.unit
    def test_affinity_with_templates(self):
        """Test affinity works alongside templates."""
        job = Job("templates_test")
        job.add_protein_chain("MVTP", ids=["A"])
        job.add_ligand(smiles="CCO", ids=["LIG"])
        
        # Add template
        job.add_template(cif="test.cif", chain_id=["A"])
        
        # Request affinity
        job.request_affinity("LIG")
        
        # Verify all components are present
        assert len(job.templates) == 1
        assert len(job.properties) == 1
        assert job.properties[0].binder == "LIG"

    @pytest.mark.unit
    def test_job_str_with_affinity(self):
        """Test that Job.__str__ correctly displays affinity properties."""
        job = Job("str_test")
        job.add_protein_chain("MVTP", ids=["A"])
        job.add_ligand(smiles="CCO", ids=["LIG"])
        job.request_affinity("LIG")
        
        job_str = str(job)
        
        assert "Properties:" in job_str
        assert "Affinity for ligand: LIG" in job_str

    @pytest.mark.unit
    def test_complete_workflow_with_warnings(self):
        """Test complete workflow with various warning conditions."""
        # Test with mixed targets to ensure warnings work in full workflow
        job = Job("complete_test")
        job.add_protein_chain("MVTPKVLH", ids=["PROT"])
        job.add_dna_chain("ATCGATCG", ids=["DNA"])
        job.add_ligand(smiles="c1ccccc1", ids=["BENZENE"])
        
        # Add various components
        job.add_pocket("BENZENE")
        job.add_template(cif="template.cif", chain_id=["PROT"])
        
        # Request affinity with expected warning
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            job.request_affinity("BENZENE")
            
            # Should warn about DNA target
            assert len(w) == 1
            assert "DnaChain" in str(w[0].message)
        
        # Verify complete job structure
        assert len(job.sequences) == 3
        assert len(job.constraints) == 1
        assert len(job.templates) == 1
        assert len(job.properties) == 1
        
        # Verify YAML generation works
        yaml_str = job.to_yaml()
        assert "affinity" in yaml_str
        assert "benzene" in yaml_str.lower() or "BENZENE" in yaml_str
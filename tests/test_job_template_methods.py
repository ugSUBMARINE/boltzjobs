"""Tests for Job class template methods (Section 2 of TODO.md)."""

import pytest
import yaml

from boltzjobs import Job
from boltzjobs.components import Template
from boltzjobs.utils import IndentedDumper


class TestJobAddTemplate:
    """Test Job.add_template() method with new parameters."""

    @pytest.fixture
    def job(self):
        """Create a basic job for testing."""
        return Job("test_job")

    @pytest.mark.unit
    def test_add_template_cif_basic(self, job):
        """Test add_template with CIF file (legacy behavior)."""
        template = job.add_template(cif="test.cif")

        assert len(job.templates) == 1
        assert isinstance(template, Template)
        assert template.cif == "test.cif"
        assert template.pdb is None
        assert template.force is False
        assert template.threshold is None

    @pytest.mark.unit
    def test_add_template_pdb_basic(self, job):
        """Test add_template with PDB file."""
        template = job.add_template(pdb="test.pdb")

        assert len(job.templates) == 1
        assert isinstance(template, Template)
        assert template.pdb == "test.pdb"
        assert template.cif is None
        assert template.force is False
        assert template.threshold is None

    @pytest.mark.unit
    def test_add_template_with_chain_ids(self, job):
        """Test add_template with chain IDs."""
        # Single chain ID as string
        template1 = job.add_template(cif="test1.cif", chain_id="A")
        assert template1.chain_id == ["A"]

        # Multiple chain IDs as list
        template2 = job.add_template(pdb="test2.pdb", chain_id=["B", "C"])
        assert template2.chain_id == ["B", "C"]

    @pytest.mark.unit
    def test_add_template_with_template_ids(self, job):
        """Test add_template with template IDs."""
        # Single template ID as string
        template1 = job.add_template(cif="test1.cif", template_id="X")
        assert template1.template_id == ["X"]

        # Multiple template IDs as list
        template2 = job.add_template(pdb="test2.pdb", template_id=["Y", "Z"])
        assert template2.template_id == ["Y", "Z"]

    @pytest.mark.unit
    def test_add_template_with_force_and_threshold(self, job):
        """Test add_template with force and threshold."""
        template = job.add_template(cif="test.cif", force=True, threshold=2.5)
        assert template.force is True
        assert template.threshold == 2.5

    @pytest.mark.unit
    def test_add_template_complete_parameters(self, job):
        """Test add_template with all parameters."""
        template = job.add_template(
            pdb="template.pdb",
            chain_id=["A", "B"],
            template_id=["X", "Y"],
            force=True,
            threshold=3.0,
        )

        assert template.pdb == "template.pdb"
        assert template.cif is None
        assert template.chain_id == ["A", "B"]
        assert template.template_id == ["X", "Y"]
        assert template.force is True
        assert template.threshold == 3.0

    @pytest.mark.unit
    def test_add_template_no_file_error(self, job):
        """Test add_template fails when no file is provided."""
        with pytest.raises(
            ValueError, match="Either 'cif' or 'pdb' file path must be provided"
        ):
            job.add_template()

    @pytest.mark.unit
    def test_add_template_both_files_error(self, job):
        """Test add_template fails when both files are provided."""
        with pytest.raises(ValueError, match="'cif' and 'pdb' are mutually exclusive"):
            job.add_template(cif="test.cif", pdb="test.pdb")

    @pytest.mark.unit
    def test_add_template_multiple_templates(self, job):
        """Test adding multiple templates to job."""
        template1 = job.add_template(cif="test1.cif")
        template2 = job.add_template(pdb="test2.pdb")
        template3 = job.add_template(cif="test3.cif", force=True, threshold=2.0)

        assert len(job.templates) == 3
        assert job.templates[0] == template1
        assert job.templates[1] == template2
        assert job.templates[2] == template3


class TestJobAddPdbTemplate:
    """Test Job.add_pdb_template() convenience method."""

    @pytest.fixture
    def job(self):
        """Create a basic job for testing."""
        return Job("test_job")

    @pytest.mark.unit
    def test_add_pdb_template_basic(self, job):
        """Test add_pdb_template basic functionality."""
        template = job.add_pdb_template("test.pdb")

        assert len(job.templates) == 1
        assert isinstance(template, Template)
        assert template.pdb == "test.pdb"
        assert template.cif is None
        assert template.force is False
        assert template.threshold is None

    @pytest.mark.unit
    def test_add_pdb_template_with_chain_ids(self, job):
        """Test add_pdb_template with chain IDs."""
        template = job.add_pdb_template("test.pdb", chain_id=["A", "B"])
        assert template.chain_id == ["A", "B"]

    @pytest.mark.unit
    def test_add_pdb_template_with_template_ids(self, job):
        """Test add_pdb_template with template IDs."""
        template = job.add_pdb_template("test.pdb", template_id=["X", "Y"])
        assert template.template_id == ["X", "Y"]

    @pytest.mark.unit
    def test_add_pdb_template_with_force(self, job):
        """Test add_pdb_template with force and threshold."""
        template = job.add_pdb_template("test.pdb", force=True, threshold=3.5)
        assert template.force is True
        assert template.threshold == 3.5

    @pytest.mark.unit
    def test_add_pdb_template_complete(self, job):
        """Test add_pdb_template with all parameters."""
        template = job.add_pdb_template(
            "template.pdb", chain_id=["A"], template_id=["X"], force=True, threshold=2.0
        )

        assert template.pdb == "template.pdb"
        assert template.chain_id == ["A"]
        assert template.template_id == ["X"]
        assert template.force is True
        assert template.threshold == 2.0


class TestJobAddCifTemplate:
    """Test Job.add_cif_template() convenience method."""

    @pytest.fixture
    def job(self):
        """Create a basic job for testing."""
        return Job("test_job")

    @pytest.mark.unit
    def test_add_cif_template_basic(self, job):
        """Test add_cif_template basic functionality."""
        template = job.add_cif_template("test.cif")

        assert len(job.templates) == 1
        assert isinstance(template, Template)
        assert template.cif == "test.cif"
        assert template.pdb is None
        assert template.force is False
        assert template.threshold is None

    @pytest.mark.unit
    def test_add_cif_template_with_chain_ids(self, job):
        """Test add_cif_template with chain IDs."""
        template = job.add_cif_template("test.cif", chain_id=["A", "B"])
        assert template.chain_id == ["A", "B"]

    @pytest.mark.unit
    def test_add_cif_template_with_template_ids(self, job):
        """Test add_cif_template with template IDs."""
        template = job.add_cif_template("test.cif", template_id=["X", "Y"])
        assert template.template_id == ["X", "Y"]

    @pytest.mark.unit
    def test_add_cif_template_with_force(self, job):
        """Test add_cif_template with force and threshold."""
        template = job.add_cif_template("test.cif", force=True, threshold=3.5)
        assert template.force is True
        assert template.threshold == 3.5

    @pytest.mark.unit
    def test_add_cif_template_complete(self, job):
        """Test add_cif_template with all parameters."""
        template = job.add_cif_template(
            "template.cif", chain_id=["A"], template_id=["X"], force=True, threshold=2.0
        )

        assert template.cif == "template.cif"
        assert template.pdb is None
        assert template.chain_id == ["A"]
        assert template.template_id == ["X"]
        assert template.force is True
        assert template.threshold == 2.0


class TestJobTemplateMethods:
    """Test integration of Job template methods."""

    @pytest.fixture
    def job_with_sequence(self):
        """Create a job with a protein sequence for testing."""
        job = Job("test_job")
        job.add_protein_chain("MVTP", ids=["A"])
        return job

    @pytest.mark.unit
    def test_job_template_yaml_serialization_cif(self, job_with_sequence):
        """Test YAML serialization with CIF template."""
        job_with_sequence.add_template(cif="test.cif", chain_id=["A"])

        yaml_dict = job_with_sequence.to_dict()
        assert "templates" in yaml_dict
        assert len(yaml_dict["templates"]) == 1

        template_dict = yaml_dict["templates"][0]
        assert template_dict["cif"] == "test.cif"
        assert template_dict["chain_id"] == "A"
        assert "pdb" not in template_dict
        assert "force" not in template_dict
        assert "threshold" not in template_dict

    @pytest.mark.unit
    def test_job_template_yaml_serialization_pdb(self, job_with_sequence):
        """Test YAML serialization with PDB template."""
        job_with_sequence.add_pdb_template("test.pdb", chain_id=["A"])

        yaml_dict = job_with_sequence.to_dict()
        template_dict = yaml_dict["templates"][0]
        assert template_dict["pdb"] == "test.pdb"
        assert template_dict["chain_id"] == "A"
        assert "cif" not in template_dict

    @pytest.mark.unit
    def test_job_template_yaml_serialization_enforced(self, job_with_sequence):
        """Test YAML serialization with enforced template."""
        job_with_sequence.add_cif_template(
            "test.cif", chain_id=["A"], force=True, threshold=2.5
        )

        yaml_dict = job_with_sequence.to_dict()
        template_dict = yaml_dict["templates"][0]
        assert template_dict["cif"] == "test.cif"
        assert template_dict["chain_id"] == "A"
        assert template_dict["force"] is True
        assert template_dict["threshold"] == 2.5

    @pytest.mark.unit
    def test_job_template_yaml_output(self, job_with_sequence):
        """Test actual YAML output format."""
        job_with_sequence.add_template(
            pdb="template.pdb",
            chain_id=["A"],
            template_id=["X"],
            force=True,
            threshold=3.0,
        )

        yaml_content = job_with_sequence.to_yaml()

        # Parse the YAML to verify structure
        parsed = yaml.safe_load(yaml_content)
        template = parsed["templates"][0]

        assert template["pdb"] == "template.pdb"
        assert template["chain_id"] == "A"
        assert template["template_id"] == "X"
        assert template["force"] is True
        assert template["threshold"] == 3.0

    @pytest.mark.unit
    def test_job_multiple_template_types(self, job_with_sequence):
        """Test job with multiple template types."""
        # Add different types of templates
        job_with_sequence.add_template(cif="basic.cif")
        job_with_sequence.add_pdb_template("standard.pdb", chain_id=["A"])
        job_with_sequence.add_cif_template(
            "enforced.cif", chain_id=["A"], force=True, threshold=2.0
        )

        assert len(job_with_sequence.templates) == 3

        yaml_dict = job_with_sequence.to_dict()
        templates = yaml_dict["templates"]

        # First template: basic CIF
        assert templates[0]["cif"] == "basic.cif"
        assert "force" not in templates[0]

        # Second template: PDB with chain ID
        assert templates[1]["pdb"] == "standard.pdb"
        assert templates[1]["chain_id"] == "A"

        # Third template: enforced CIF
        assert templates[2]["cif"] == "enforced.cif"
        assert templates[2]["force"] is True
        assert templates[2]["threshold"] == 2.0


class TestJobTemplateBackwardCompatibility:
    """Test backward compatibility of Job template methods."""

    @pytest.fixture
    def job_with_sequence(self):
        """Create a job with a protein sequence for testing."""
        job = Job("test_job")
        job.add_protein_chain("MVTP", ids=["A"])
        return job

    @pytest.mark.unit
    def test_legacy_add_template_still_works(self, job_with_sequence):
        """Test that legacy add_template calls still work."""
        # This is how add_template was called before the changes
        template = job_with_sequence.add_template(
            cif="legacy.cif", chain_id="A", template_id="X"
        )

        assert template.cif == "legacy.cif"
        assert template.pdb is None
        assert template.chain_id == ["A"]
        assert template.template_id == ["X"]
        assert template.force is False
        assert template.threshold is None

    @pytest.mark.unit
    def test_legacy_template_serialization(self, job_with_sequence):
        """Test that legacy templates serialize correctly."""
        job_with_sequence.add_template(cif="legacy.cif", chain_id="A")

        yaml_dict = job_with_sequence.to_dict()
        template_dict = yaml_dict["templates"][0]

        # Should match old format
        assert template_dict == {"cif": "legacy.cif", "chain_id": "A"}

    @pytest.mark.unit
    def test_yaml_round_trip_legacy(self, job_with_sequence):
        """Test round-trip YAML serialization for legacy templates."""
        job_with_sequence.add_template(cif="legacy.cif", chain_id=["A", "B"])

        # Serialize to YAML
        yaml_content = job_with_sequence.to_yaml()

        # Parse back from YAML
        yaml_dict = yaml.safe_load(yaml_content)

        template_data = yaml_dict["templates"][0]
        assert template_data["cif"] == "legacy.cif"
        assert template_data["chain_id"] == ["A", "B"]
        assert "pdb" not in template_data
        assert "force" not in template_data
        assert "threshold" not in template_data


class TestJobTemplateIntegration:
    """Test template integration with full job workflow."""

    @pytest.mark.unit
    def test_job_with_templates_complete_workflow(self):
        """Test complete job workflow with templates."""
        job = Job("complex_job")

        # Add sequences
        job.add_protein_chain("MVTPKVLH", ids=["A"])
        job.add_ligand(smiles="CCO", ids=["LIG"])

        # Add various templates
        job.add_template(cif="protein_template.cif", chain_id=["A"])
        job.add_pdb_template("backup_template.pdb", chain_id=["A"])
        job.add_cif_template(
            "critical_template.cif",
            chain_id=["A"],
            template_id=["P"],
            force=True,
            threshold=2.0,
        )

        # Add constraints and affinity
        job.add_pocket("LIG")
        job.request_affinity("LIG")

        # Verify job structure
        assert len(job.sequences) == 2
        assert len(job.templates) == 3
        assert len(job.constraints) == 1
        assert len(job.properties) == 1

        # Test serialization
        yaml_dict = job.to_dict()
        assert "sequences" in yaml_dict
        assert "templates" in yaml_dict
        assert "constraints" in yaml_dict
        assert "properties" in yaml_dict

        # Verify template structure
        templates = yaml_dict["templates"]
        assert len(templates) == 3

        assert templates[0]["cif"] == "protein_template.cif"
        assert templates[1]["pdb"] == "backup_template.pdb"
        assert templates[2]["cif"] == "critical_template.cif"
        assert templates[2]["force"] is True
        assert templates[2]["threshold"] == 2.0

    @pytest.mark.unit
    def test_job_template_yaml_output_format(self):
        """Test that template YAML output matches expected schema."""
        job = Job("schema_test")
        job.add_protein_chain("MVTP", ids=["A"])

        # Add template with all features
        job.add_template(
            pdb="test.pdb",
            chain_id=["A", "B"],
            template_id=["P1", "P2"],
            force=True,
            threshold=2.5,
        )

        yaml_content = job.to_yaml()

        # Verify YAML structure matches expected Boltz-2 schema
        expected_lines = [
            "templates:",
            "- pdb: test.pdb",
            "  chain_id: [A, B]",
            "  template_id: [P1, P2]",
            "  force: true",
            "  threshold: 2.5",
        ]

        for expected_line in expected_lines:
            assert expected_line in yaml_content

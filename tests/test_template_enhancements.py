"""Tests for Template class enhancements (Section 2 of TODO.md)."""

import pytest
from dataclasses import FrozenInstanceError

from boltzjobs.components import Template
from boltzjobs.utils import FlowStyleList


class TestTemplateEnhancements:
    """Test Template class with PDB support, force and threshold parameters."""

    @pytest.mark.unit
    def test_template_cif_creation(self):
        """Test Template creation with CIF file."""
        template = Template(cif="test.cif")
        assert template.cif == "test.cif"
        assert template.pdb is None
        assert template.force is False
        assert template.threshold is None
        assert template.chain_id == []
        assert template.template_id == []

    @pytest.mark.unit
    def test_template_pdb_creation(self):
        """Test Template creation with PDB file."""
        template = Template(pdb="test.pdb")
        assert template.pdb == "test.pdb"
        assert template.cif is None
        assert template.force is False
        assert template.threshold is None
        assert template.chain_id == []
        assert template.template_id == []

    @pytest.mark.unit
    def test_template_with_all_parameters(self):
        """Test Template creation with all parameters."""
        template = Template(
            pdb="template.pdb",
            chain_id=["A", "B"],
            template_id=["X", "Y"],
            force=True,
            threshold=2.5
        )
        assert template.pdb == "template.pdb"
        assert template.cif is None
        assert template.chain_id == ["A", "B"]
        assert template.template_id == ["X", "Y"]
        assert template.force is True
        assert template.threshold == 2.5

    @pytest.mark.unit
    def test_template_no_file_error(self):
        """Test Template creation fails when no file is provided."""
        with pytest.raises(ValueError, match="Either 'cif' or 'pdb' file path must be provided"):
            Template()

    @pytest.mark.unit
    def test_template_both_files_error(self):
        """Test Template creation fails when both files are provided."""
        with pytest.raises(ValueError, match="'cif' and 'pdb' are mutually exclusive"):
            Template(cif="test.cif", pdb="test.pdb")

    @pytest.mark.unit
    def test_template_force_without_threshold_error(self):
        """Test Template creation fails when force=True but no threshold."""
        with pytest.raises(ValueError, match="'threshold' must be specified when 'force=True'"):
            Template(cif="test.cif", force=True)

    @pytest.mark.unit
    def test_template_force_with_threshold_valid(self):
        """Test Template creation succeeds when force=True and threshold provided."""
        template = Template(cif="test.cif", force=True, threshold=3.0)
        assert template.force is True
        assert template.threshold == 3.0

    @pytest.mark.unit
    def test_template_threshold_without_force_valid(self):
        """Test Template creation allows threshold without force=True."""
        template = Template(cif="test.cif", threshold=2.0)
        assert template.force is False
        assert template.threshold == 2.0


class TestTemplateStringRepresentation:
    """Test Template string representation with new fields."""

    @pytest.mark.unit
    def test_template_str_cif_basic(self):
        """Test string representation for CIF template."""
        template = Template(cif="test.cif")
        str_repr = str(template)
        assert "++ Template: CIF: 'test.cif'" in str_repr

    @pytest.mark.unit
    def test_template_str_pdb_basic(self):
        """Test string representation for PDB template."""
        template = Template(pdb="test.pdb")
        str_repr = str(template)
        assert "++ Template: PDB: 'test.pdb'" in str_repr

    @pytest.mark.unit
    def test_template_str_with_chain_ids(self):
        """Test string representation with chain IDs."""
        template = Template(cif="test.cif", chain_id=["A", "B"])
        str_repr = str(template)
        assert "For chain ID(s): A, B" in str_repr

    @pytest.mark.unit
    def test_template_str_with_template_ids(self):
        """Test string representation with template IDs."""
        template = Template(pdb="test.pdb", template_id=["X", "Y"])
        str_repr = str(template)
        assert "Template ID(s): X, Y" in str_repr

    @pytest.mark.unit
    def test_template_str_with_force_and_threshold(self):
        """Test string representation with force and threshold."""
        template = Template(cif="test.cif", force=True, threshold=2.5)
        str_repr = str(template)
        assert "Force: True (threshold: 2.5Å)" in str_repr

    @pytest.mark.unit
    def test_template_str_complete(self):
        """Test string representation with all fields."""
        template = Template(
            pdb="template.pdb",
            chain_id=["A"],
            template_id=["X"],
            force=True,
            threshold=3.0
        )
        str_repr = str(template)
        assert "++ Template: PDB: 'template.pdb'" in str_repr
        assert "Template ID(s): X" in str_repr
        assert "For chain ID(s): A" in str_repr
        assert "Force: True (threshold: 3.0Å)" in str_repr


class TestTemplateSerialization:
    """Test Template serialization and deserialization."""

    @pytest.mark.unit
    def test_template_to_dict_cif_basic(self):
        """Test Template.to_dict() for CIF file."""
        template = Template(cif="test.cif")
        result = template.to_dict()
        assert result == {"cif": "test.cif"}

    @pytest.mark.unit
    def test_template_to_dict_pdb_basic(self):
        """Test Template.to_dict() for PDB file."""
        template = Template(pdb="test.pdb")
        result = template.to_dict()
        assert result == {"pdb": "test.pdb"}

    @pytest.mark.unit
    def test_template_to_dict_with_single_chain_id(self):
        """Test Template.to_dict() with single chain ID."""
        template = Template(cif="test.cif", chain_id=["A"])
        result = template.to_dict()
        assert result == {"cif": "test.cif", "chain_id": "A"}

    @pytest.mark.unit
    def test_template_to_dict_with_multiple_chain_ids(self):
        """Test Template.to_dict() with multiple chain IDs."""
        template = Template(cif="test.cif", chain_id=["A", "B"])
        result = template.to_dict()
        assert result["cif"] == "test.cif"
        assert isinstance(result["chain_id"], FlowStyleList)
        assert result["chain_id"] == ["A", "B"]

    @pytest.mark.unit
    def test_template_to_dict_with_single_template_id(self):
        """Test Template.to_dict() with single template ID."""
        template = Template(pdb="test.pdb", template_id=["X"])
        result = template.to_dict()
        assert result == {"pdb": "test.pdb", "template_id": "X"}

    @pytest.mark.unit
    def test_template_to_dict_with_multiple_template_ids(self):
        """Test Template.to_dict() with multiple template IDs."""
        template = Template(pdb="test.pdb", template_id=["X", "Y"])
        result = template.to_dict()
        assert result["pdb"] == "test.pdb"
        assert isinstance(result["template_id"], FlowStyleList)
        assert result["template_id"] == ["X", "Y"]

    @pytest.mark.unit
    def test_template_to_dict_with_force_true(self):
        """Test Template.to_dict() with force=True."""
        template = Template(cif="test.cif", force=True, threshold=2.5)
        result = template.to_dict()
        expected = {"cif": "test.cif", "force": True, "threshold": 2.5}
        assert result == expected

    @pytest.mark.unit
    def test_template_to_dict_with_force_false(self):
        """Test Template.to_dict() with force=False (default)."""
        template = Template(cif="test.cif")
        result = template.to_dict()
        # force=False should not appear in dict (default value)
        assert "force" not in result
        assert result == {"cif": "test.cif"}

    @pytest.mark.unit
    def test_template_to_dict_with_threshold_only(self):
        """Test Template.to_dict() with threshold but no force."""
        template = Template(cif="test.cif", threshold=3.0)
        result = template.to_dict()
        expected = {"cif": "test.cif", "threshold": 3.0}
        assert result == expected

    @pytest.mark.unit
    def test_template_to_dict_complete(self):
        """Test Template.to_dict() with all fields."""
        template = Template(
            pdb="template.pdb",
            chain_id=["A", "B"],
            template_id=["X", "Y"],
            force=True,
            threshold=2.0
        )
        result = template.to_dict()
        
        expected = {
            "pdb": "template.pdb",
            "chain_id": FlowStyleList(["A", "B"]),
            "template_id": FlowStyleList(["X", "Y"]),
            "force": True,
            "threshold": 2.0
        }
        assert result == expected


class TestTemplateDeserialization:
    """Test Template.from_dict() method."""

    @pytest.mark.unit
    def test_template_from_dict_cif_basic(self):
        """Test Template.from_dict() for CIF file."""
        data = {"cif": "test.cif"}
        template = Template.from_dict(data)
        assert template.cif == "test.cif"
        assert template.pdb is None
        assert template.force is False
        assert template.threshold is None

    @pytest.mark.unit
    def test_template_from_dict_pdb_basic(self):
        """Test Template.from_dict() for PDB file."""
        data = {"pdb": "test.pdb"}
        template = Template.from_dict(data)
        assert template.pdb == "test.pdb"
        assert template.cif is None
        assert template.force is False
        assert template.threshold is None

    @pytest.mark.unit
    def test_template_from_dict_with_single_chain_id(self):
        """Test Template.from_dict() with single chain ID."""
        data = {"cif": "test.cif", "chain_id": "A"}
        template = Template.from_dict(data)
        assert template.chain_id == ["A"]

    @pytest.mark.unit
    def test_template_from_dict_with_multiple_chain_ids(self):
        """Test Template.from_dict() with multiple chain IDs."""
        data = {"cif": "test.cif", "chain_id": ["A", "B"]}
        template = Template.from_dict(data)
        assert template.chain_id == ["A", "B"]

    @pytest.mark.unit
    def test_template_from_dict_with_single_template_id(self):
        """Test Template.from_dict() with single template ID."""
        data = {"pdb": "test.pdb", "template_id": "X"}
        template = Template.from_dict(data)
        assert template.template_id == ["X"]

    @pytest.mark.unit
    def test_template_from_dict_with_multiple_template_ids(self):
        """Test Template.from_dict() with multiple template IDs."""
        data = {"pdb": "test.pdb", "template_id": ["X", "Y"]}
        template = Template.from_dict(data)
        assert template.template_id == ["X", "Y"]

    @pytest.mark.unit
    def test_template_from_dict_with_force_and_threshold(self):
        """Test Template.from_dict() with force and threshold."""
        data = {"cif": "test.cif", "force": True, "threshold": 2.5}
        template = Template.from_dict(data)
        assert template.force is True
        assert template.threshold == 2.5

    @pytest.mark.unit
    def test_template_from_dict_with_defaults(self):
        """Test Template.from_dict() uses proper defaults."""
        data = {"pdb": "test.pdb"}
        template = Template.from_dict(data)
        assert template.chain_id == []
        assert template.template_id == []
        assert template.force is False
        assert template.threshold is None

    @pytest.mark.unit
    def test_template_from_dict_complete(self):
        """Test Template.from_dict() with all fields."""
        data = {
            "pdb": "template.pdb",
            "chain_id": ["A", "B"],
            "template_id": ["X", "Y"],
            "force": True,
            "threshold": 2.0
        }
        template = Template.from_dict(data)
        assert template.pdb == "template.pdb"
        assert template.cif is None
        assert template.chain_id == ["A", "B"]
        assert template.template_id == ["X", "Y"]
        assert template.force is True
        assert template.threshold == 2.0


class TestTemplateRoundTripSerialization:
    """Test Template round-trip serialization."""

    @pytest.mark.unit
    def test_template_round_trip_cif_basic(self):
        """Test Template round-trip serialization for CIF."""
        original = Template(cif="test.cif")
        data = original.to_dict()
        restored = Template.from_dict(data)
        assert restored.cif == original.cif
        assert restored.pdb == original.pdb
        assert restored.force == original.force
        assert restored.threshold == original.threshold

    @pytest.mark.unit
    def test_template_round_trip_pdb_basic(self):
        """Test Template round-trip serialization for PDB."""
        original = Template(pdb="test.pdb")
        data = original.to_dict()
        restored = Template.from_dict(data)
        assert restored.pdb == original.pdb
        assert restored.cif == original.cif
        assert restored.force == original.force
        assert restored.threshold == original.threshold

    @pytest.mark.unit
    def test_template_round_trip_complete(self):
        """Test Template round-trip serialization with all fields."""
        original = Template(
            pdb="template.pdb",
            chain_id=["A", "B"],
            template_id=["X", "Y"],
            force=True,
            threshold=2.5
        )
        data = original.to_dict()
        restored = Template.from_dict(data)
        
        assert restored.pdb == original.pdb
        assert restored.cif == original.cif
        assert restored.chain_id == original.chain_id
        assert restored.template_id == original.template_id
        assert restored.force == original.force
        assert restored.threshold == original.threshold


class TestTemplateBackwardCompatibility:
    """Test backward compatibility with existing Template usage."""

    @pytest.mark.unit
    def test_template_legacy_cif_creation(self):
        """Test Template still works with legacy CIF-only creation."""
        # This mimics the old Template(cif, chain_id=..., template_id=...)
        template = Template(cif="legacy.cif", chain_id=["A"], template_id=["X"])
        assert template.cif == "legacy.cif"
        assert template.pdb is None
        assert template.chain_id == ["A"]
        assert template.template_id == ["X"]
        assert template.force is False
        assert template.threshold is None

    @pytest.mark.unit
    def test_template_legacy_serialization(self):
        """Test Template serialization is backward compatible."""
        template = Template(cif="legacy.cif", chain_id=["A"])
        result = template.to_dict()
        
        # Should look like old format but with new structure
        expected = {"cif": "legacy.cif", "chain_id": "A"}
        assert result == expected

    @pytest.mark.unit
    def test_template_legacy_deserialization(self):
        """Test Template can deserialize old-format data."""
        # Old format data
        data = {"cif": "legacy.cif", "chain_id": "A", "template_id": ["X", "Y"]}
        template = Template.from_dict(data)
        
        assert template.cif == "legacy.cif"
        assert template.pdb is None
        assert template.chain_id == ["A"]
        assert template.template_id == ["X", "Y"]
        assert template.force is False
        assert template.threshold is None


class TestTemplateEdgeCases:
    """Test Template edge cases and error conditions."""

    @pytest.mark.unit
    def test_template_empty_string_cif(self):
        """Test Template with empty string CIF is invalid."""
        with pytest.raises(ValueError, match="Either 'cif' or 'pdb' file path must be provided"):
            Template(cif="")

    @pytest.mark.unit
    def test_template_empty_string_pdb(self):
        """Test Template with empty string PDB is invalid."""
        with pytest.raises(ValueError, match="Either 'cif' or 'pdb' file path must be provided"):
            Template(pdb="")

    @pytest.mark.unit
    def test_template_none_values_both(self):
        """Test Template with explicit None values for both files."""
        with pytest.raises(ValueError, match="Either 'cif' or 'pdb' file path must be provided"):
            Template(cif=None, pdb=None)

    @pytest.mark.unit
    def test_template_empty_strings_both(self):
        """Test Template with empty strings for both files."""
        # Empty strings are treated as falsy, so this is equivalent to no files provided
        with pytest.raises(ValueError, match="Either 'cif' or 'pdb' file path must be provided"):
            Template(cif="", pdb="")

    @pytest.mark.unit
    def test_template_force_true_threshold_zero(self):
        """Test Template allows threshold=0.0 with force=True."""
        template = Template(cif="test.cif", force=True, threshold=0.0)
        assert template.threshold == 0.0

    @pytest.mark.unit
    def test_template_force_true_threshold_negative(self):
        """Test Template allows negative threshold with force=True."""
        template = Template(cif="test.cif", force=True, threshold=-1.0)
        assert template.threshold == -1.0
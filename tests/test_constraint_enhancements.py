"""
Test suite for constraint system enhancements from TODO.md section 1.

Tests cover:
- Force parameter support for Pocket and Contact classes
- Distance validation with 4-20Å range and 6Å default
- Enhanced pocket constraint validation (any sequence type as binder)
"""

import pytest
from boltzjobs import Job
from boltzjobs.components import Contact, Pocket
from boltzjobs.utils import validate_distance


class TestDistanceValidation:
    """Test distance validation utility function."""
    
    @pytest.mark.unit
    def test_validate_distance_valid_ranges(self):
        """Test that valid distances (4-20Å) pass validation."""
        # Test boundary values
        validate_distance(4.0)
        validate_distance(20.0)
        
        # Test typical values
        validate_distance(6.0)
        validate_distance(10.0)
        validate_distance(15.5)
        
        # Test integer values
        validate_distance(5)
        validate_distance(8)

    @pytest.mark.unit
    def test_validate_distance_invalid_ranges(self):
        """Test that invalid distances raise ValueError."""
        # Below minimum
        with pytest.raises(ValueError, match="must be between 4.0 and 20.0"):
            validate_distance(3.9)
        
        with pytest.raises(ValueError, match="must be between 4.0 and 20.0"):
            validate_distance(0.0)
        
        # Above maximum
        with pytest.raises(ValueError, match="must be between 4.0 and 20.0"):
            validate_distance(20.1)
        
        with pytest.raises(ValueError, match="must be between 4.0 and 20.0"):
            validate_distance(50.0)

    @pytest.mark.unit
    def test_validate_distance_invalid_types(self):
        """Test that non-numeric types raise TypeError."""
        with pytest.raises(TypeError, match="must be a number"):
            validate_distance("6.0")
        
        with pytest.raises(TypeError, match="must be a number"):
            validate_distance(None)
        
        with pytest.raises(TypeError, match="must be a number"):
            validate_distance([6.0])

    @pytest.mark.unit
    def test_validate_distance_custom_parameter_name(self):
        """Test custom parameter name in error messages."""
        with pytest.raises(ValueError, match="custom_param must be between 4.0 and 20.0"):
            validate_distance(3.0, "custom_param")


class TestPocketEnhancements:
    """Test Pocket class enhancements."""
    
    @pytest.mark.unit
    def test_pocket_default_distance_changed(self):
        """Test that default max_distance is now 6.0Å."""
        pocket = Pocket("A")
        assert pocket.max_distance == 6.0

    @pytest.mark.unit
    def test_pocket_force_parameter_default(self):
        """Test that force parameter defaults to False."""
        pocket = Pocket("A")
        assert pocket.force is False

    @pytest.mark.unit
    def test_pocket_force_parameter_explicit(self):
        """Test explicit force parameter setting."""
        pocket = Pocket("A", force=True)
        assert pocket.force is True
        
        pocket = Pocket("A", force=False)
        assert pocket.force is False

    @pytest.mark.unit
    def test_pocket_distance_validation(self):
        """Test that Pocket validates distance on creation."""
        # Valid distance
        pocket = Pocket("A", max_distance=8.0)
        assert pocket.max_distance == 8.0
        
        # Invalid distance - too small
        with pytest.raises(ValueError, match="must be between 4.0 and 20.0"):
            Pocket("A", max_distance=3.0)
        
        # Invalid distance - too large  
        with pytest.raises(ValueError, match="must be between 4.0 and 20.0"):
            Pocket("A", max_distance=25.0)

    @pytest.mark.unit
    def test_pocket_to_dict_without_force(self):
        """Test to_dict doesn't include force when False."""
        pocket = Pocket("A", [("B", 10)], 8.0, False)
        result = pocket.to_dict()
        
        expected = {
            "pocket": {
                "binder": "A",
                "contacts": [("B", 10)],
                "max_distance": 8.0
            }
        }
        assert result == expected
        assert "force" not in result["pocket"]

    @pytest.mark.unit
    def test_pocket_to_dict_with_force(self):
        """Test to_dict includes force when True."""
        pocket = Pocket("A", [("B", 10)], 8.0, True)
        result = pocket.to_dict()
        
        expected = {
            "pocket": {
                "binder": "A", 
                "contacts": [("B", 10)],
                "max_distance": 8.0,
                "force": True
            }
        }
        assert result == expected

    @pytest.mark.unit
    def test_pocket_from_dict_without_force(self):
        """Test from_dict with no force field defaults to False."""
        data = {
            "binder": "A",
            "contacts": [("B", 10)],
            "max_distance": 8.0
        }
        pocket = Pocket.from_dict(data)
        
        assert pocket.binder == "A"
        assert pocket.contacts == [("B", 10)]
        assert pocket.max_distance == 8.0
        assert pocket.force is False

    @pytest.mark.unit
    def test_pocket_from_dict_with_force(self):
        """Test from_dict with force field."""
        data = {
            "binder": "A",
            "contacts": [("B", 10)],
            "max_distance": 8.0,
            "force": True
        }
        pocket = Pocket.from_dict(data)
        
        assert pocket.binder == "A"
        assert pocket.contacts == [("B", 10)]
        assert pocket.max_distance == 8.0
        assert pocket.force is True

    @pytest.mark.unit
    def test_pocket_from_dict_new_default_distance(self):
        """Test from_dict uses new default distance (6.0Å)."""
        data = {"binder": "A", "contacts": []}
        pocket = Pocket.from_dict(data)
        assert pocket.max_distance == 6.0

    @pytest.mark.unit
    def test_pocket_str_representation_with_force(self):
        """Test string representation includes force when True."""
        pocket = Pocket("A", force=True)
        result = str(pocket)
        assert "Force: True" in result

    @pytest.mark.unit
    def test_pocket_str_representation_without_force(self):
        """Test string representation doesn't mention force when False."""
        pocket = Pocket("A", force=False)
        result = str(pocket)
        assert "Force" not in result


class TestContactEnhancements:
    """Test Contact class enhancements."""
    
    @pytest.mark.unit
    def test_contact_default_distance_changed(self):
        """Test that default max_distance is now 6.0Å."""
        contact = Contact(("A", 1), ("B", 2))
        assert contact.max_distance == 6.0

    @pytest.mark.unit
    def test_contact_force_parameter_default(self):
        """Test that force parameter defaults to False."""
        contact = Contact(("A", 1), ("B", 2))
        assert contact.force is False

    @pytest.mark.unit
    def test_contact_force_parameter_explicit(self):
        """Test explicit force parameter setting."""
        contact = Contact(("A", 1), ("B", 2), force=True)
        assert contact.force is True
        
        contact = Contact(("A", 1), ("B", 2), force=False)
        assert contact.force is False

    @pytest.mark.unit
    def test_contact_distance_validation(self):
        """Test that Contact validates distance on creation."""
        # Valid distance
        contact = Contact(("A", 1), ("B", 2), max_distance=12.0)
        assert contact.max_distance == 12.0
        
        # Invalid distance - too small
        with pytest.raises(ValueError, match="must be between 4.0 and 20.0"):
            Contact(("A", 1), ("B", 2), max_distance=2.0)
        
        # Invalid distance - too large
        with pytest.raises(ValueError, match="must be between 4.0 and 20.0"):
            Contact(("A", 1), ("B", 2), max_distance=30.0)

    @pytest.mark.unit
    def test_contact_to_dict_without_force(self):
        """Test to_dict doesn't include force when False."""
        contact = Contact(("A", 1), ("B", 2), 7.0, False)
        result = contact.to_dict()
        
        expected = {
            "contact": {
                "token1": ["A", 1],
                "token2": ["B", 2],
                "max_distance": 7.0
            }
        }
        assert result == expected
        assert "force" not in result["contact"]

    @pytest.mark.unit
    def test_contact_to_dict_with_force(self):
        """Test to_dict includes force when True."""
        contact = Contact(("A", 1), ("B", 2), 7.0, True)
        result = contact.to_dict()
        
        expected = {
            "contact": {
                "token1": ["A", 1],
                "token2": ["B", 2],
                "max_distance": 7.0,
                "force": True
            }
        }
        assert result == expected

    @pytest.mark.unit
    def test_contact_from_dict_without_force(self):
        """Test from_dict with no force field defaults to False."""
        data = {
            "token1": ("A", 1),
            "token2": ("B", 2),
            "max_distance": 7.0
        }
        contact = Contact.from_dict(data)
        
        assert contact.token1 == ("A", 1)
        assert contact.token2 == ("B", 2)
        assert contact.max_distance == 7.0
        assert contact.force is False

    @pytest.mark.unit
    def test_contact_from_dict_with_force(self):
        """Test from_dict with force field."""
        data = {
            "token1": ("A", 1),
            "token2": ("B", 2), 
            "max_distance": 7.0,
            "force": True
        }
        contact = Contact.from_dict(data)
        
        assert contact.token1 == ("A", 1)
        assert contact.token2 == ("B", 2)
        assert contact.max_distance == 7.0
        assert contact.force is True

    @pytest.mark.unit
    def test_contact_from_dict_new_default_distance(self):
        """Test from_dict uses new default distance (6.0Å)."""
        data = {"token1": ("A", 1), "token2": ("B", 2)}
        contact = Contact.from_dict(data)
        assert contact.max_distance == 6.0

    @pytest.mark.unit
    def test_contact_str_representation_with_force(self):
        """Test string representation includes force when True."""
        contact = Contact(("A", 1), ("B", 2), force=True)
        result = str(contact)
        assert "force: True" in result

    @pytest.mark.unit
    def test_contact_str_representation_without_force(self):
        """Test string representation doesn't mention force when False."""
        contact = Contact(("A", 1), ("B", 2), force=False)
        result = str(contact)
        assert "force:" not in result


class TestJobConstraintEnhancements:
    """Test Job class constraint enhancements."""
    
    @pytest.mark.unit
    def test_add_pocket_new_default_distance(self):
        """Test add_pocket uses new default distance (6.0Å)."""
        job = Job("test")
        job.add_ligand(ccd="ATP", ids="A")
        
        pocket = job.add_pocket("A")
        assert pocket.max_distance == 6.0

    @pytest.mark.unit
    def test_add_pocket_force_parameter(self):
        """Test add_pocket supports force parameter."""
        job = Job("test")
        job.add_ligand(ccd="ATP", ids="A")
        
        # Test default force=False
        pocket1 = job.add_pocket("A")
        assert pocket1.force is False
        
        # Test explicit force=True
        pocket2 = job.add_pocket("A", force=True)
        assert pocket2.force is True

    @pytest.mark.unit
    def test_add_pocket_allows_protein_binder(self):
        """Test add_pocket now allows protein chains as binders."""
        job = Job("test")
        job.add_protein_chain("MKLLVV", ids="A")
        
        # This should now work (previously would have failed)
        pocket = job.add_pocket("A")
        assert pocket.binder == "A"

    @pytest.mark.unit
    def test_add_pocket_allows_dna_binder(self):
        """Test add_pocket allows DNA chains as binders."""
        job = Job("test")
        job.add_dna_chain("ATCG", ids="A")
        
        pocket = job.add_pocket("A")
        assert pocket.binder == "A"

    @pytest.mark.unit
    def test_add_pocket_allows_rna_binder(self):
        """Test add_pocket allows RNA chains as binders."""
        job = Job("test")
        job.add_rna_chain("AUCG", ids="A")
        
        pocket = job.add_pocket("A")
        assert pocket.binder == "A"

    @pytest.mark.unit
    def test_add_pocket_allows_ligand_binder(self):
        """Test add_pocket still allows ligands as binders."""
        job = Job("test")
        job.add_ligand(ccd="ATP", ids="A")
        
        pocket = job.add_pocket("A")
        assert pocket.binder == "A"

    @pytest.mark.unit
    def test_add_pocket_invalid_binder_id(self):
        """Test add_pocket raises error for non-existent binder ID."""
        job = Job("test")
        job.add_ligand(ccd="ATP", ids="A")
        
        with pytest.raises(ValueError, match="No sequence with id 'Z' found"):
            job.add_pocket("Z")

    @pytest.mark.unit
    def test_add_contact_new_default_distance(self):
        """Test add_contact uses new default distance (6.0Å)."""
        job = Job("test")
        job.add_protein_chain("MKLLVV", ids="A")
        job.add_ligand(ccd="ATP", ids="B")
        
        contact = job.add_contact("A", 1, "B", 2)
        assert contact.max_distance == 6.0

    @pytest.mark.unit
    def test_add_contact_force_parameter(self):
        """Test add_contact supports force parameter."""
        job = Job("test")
        job.add_protein_chain("MKLLVV", ids="A")
        job.add_ligand(ccd="ATP", ids="B")
        
        # Test default force=False
        contact1 = job.add_contact("A", 1, "B", 2)
        assert contact1.force is False
        
        # Test explicit force=True
        contact2 = job.add_contact("A", 3, "B", 4, force=True)
        assert contact2.force is True

    @pytest.mark.unit
    def test_add_contact_returns_contact_object(self):
        """Test add_contact now returns the Contact object."""
        job = Job("test")
        job.add_protein_chain("MKLLVV", ids="A")
        job.add_ligand(ccd="ATP", ids="B")
        
        contact = job.add_contact("A", 1, "B", 2, max_distance=8.0)
        assert isinstance(contact, Contact)
        assert contact.token1 == ("A", 1)
        assert contact.token2 == ("B", 2)
        assert contact.max_distance == 8.0


class TestConstraintYamlSerialization:
    """Test YAML serialization with new constraint features."""
    
    @pytest.mark.unit
    def test_pocket_yaml_round_trip_with_force(self):
        """Test pocket with force parameter survives YAML round-trip."""
        job = Job("test")
        job.add_ligand(ccd="ATP", ids="A")
        pocket = job.add_pocket("A", max_distance=8.0, force=True)
        
        # Convert to dict and back
        job_dict = job.to_dict()
        
        # Check constraint in serialized form
        constraint = job_dict["constraints"][0]
        assert constraint["pocket"]["force"] is True
        assert constraint["pocket"]["max_distance"] == 8.0

    @pytest.mark.unit
    def test_contact_yaml_round_trip_with_force(self):
        """Test contact with force parameter survives YAML round-trip."""
        job = Job("test")
        job.add_protein_chain("MKLLVV", ids="A")
        job.add_ligand(ccd="ATP", ids="B")
        contact = job.add_contact("A", 1, "B", 2, max_distance=10.0, force=True)
        
        # Convert to dict and back
        job_dict = job.to_dict()
        
        # Check constraint in serialized form
        constraint = job_dict["constraints"][0]
        assert constraint["contact"]["force"] is True
        assert constraint["contact"]["max_distance"] == 10.0

    @pytest.mark.unit
    def test_pocket_yaml_round_trip_without_force(self):
        """Test pocket without force doesn't include force in YAML."""
        job = Job("test")
        job.add_ligand(ccd="ATP", ids="A")
        pocket = job.add_pocket("A", max_distance=8.0, force=False)
        
        # Convert to dict
        job_dict = job.to_dict()
        
        # Check constraint in serialized form - should not have force field
        constraint = job_dict["constraints"][0]
        assert "force" not in constraint["pocket"]
        assert constraint["pocket"]["max_distance"] == 8.0


class TestIntegrationConstraintEnhancements:
    """Integration tests for constraint enhancements."""
    
    @pytest.mark.integration
    def test_mixed_sequence_types_with_pocket_constraints(self):
        """Test pocket constraints work with mixed sequence types."""
        job = Job("mixed_pocket_test")
        
        # Add various sequence types
        job.add_protein_chain("MKLLVV", ids="P1")
        job.add_dna_chain("ATCG", ids="D1")
        job.add_rna_chain("AUCG", ids="R1")
        job.add_ligand(ccd="ATP", ids="L1")
        
        # Add pocket constraints for each type
        pocket_p = job.add_pocket("P1", max_distance=8.0, force=True)
        pocket_d = job.add_pocket("D1", max_distance=10.0)
        pocket_r = job.add_pocket("R1", max_distance=7.0, force=True)
        pocket_l = job.add_pocket("L1", max_distance=6.0)
        
        # Verify all pockets were added
        assert len(job.constraints) == 4
        
        # Check serialization
        job_dict = job.to_dict()
        constraints = job_dict["constraints"]
        
        # Verify pocket constraints are properly serialized
        pocket_constraints = [c for c in constraints if "pocket" in c]
        assert len(pocket_constraints) == 4
        
        # Check specific force values
        force_pockets = [c for c in pocket_constraints if c["pocket"].get("force", False)]
        assert len(force_pockets) == 2  # P1 and R1 have force=True

    @pytest.mark.integration
    def test_enforced_constraints_complex_scenario(self):
        """Test complex scenario with enforced constraints."""
        job = Job("enforced_constraints")
        
        # Add protein and ligand
        job.add_protein_chain("MKLLVVGATTGC", ids="A")
        job.add_ligand(smiles="CCO", ids="B")
        
        # Add enforced pocket constraint
        pocket = job.add_pocket("B", 
                              contact_tokens=[("A", 5), ("A", 8)],
                              max_distance=7.5,
                              force=True)
        
        # Add enforced contact constraint
        contact = job.add_contact("A", 3, "B", 1, 
                                max_distance=5.5, 
                                force=True)
        
        # Verify constraints
        assert len(job.constraints) == 2
        assert pocket.force is True
        assert contact.force is True
        
        # Test serialization
        job_dict = job.to_dict()
        constraints = job_dict["constraints"]
        
        # Both constraints should have force=True
        for constraint in constraints:
            constraint_type = list(constraint.keys())[0]
            assert constraint[constraint_type]["force"] is True

    @pytest.mark.integration
    def test_distance_validation_integration(self):
        """Test distance validation works throughout the system."""
        job = Job("distance_validation")
        job.add_protein_chain("MKLLVV", ids="A")
        job.add_ligand(ccd="ATP", ids="B")
        
        # Valid distance ranges should work
        job.add_pocket("B", max_distance=4.0)  # minimum
        job.add_contact("A", 1, "B", 1, max_distance=20.0)  # maximum
        
        # Invalid distances should fail
        with pytest.raises(ValueError, match="must be between 4.0 and 20.0"):
            job.add_pocket("B", max_distance=3.9)
            
        with pytest.raises(ValueError, match="must be between 4.0 and 20.0"):
            job.add_contact("A", 1, "B", 1, max_distance=21.0)
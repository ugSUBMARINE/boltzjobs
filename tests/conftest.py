"""Shared pytest fixtures for boltzjobs tests."""

import pytest
import tempfile
from pathlib import Path

from boltzjobs import Job
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
)


@pytest.fixture
def sample_protein_sequence():
    """Sample protein sequence for testing."""
    return "MVTPEGNVSLVDESLLVGVTDED"


@pytest.fixture
def sample_dna_sequence():
    """Sample DNA sequence for testing."""
    return "ATCGATCGATCG"


@pytest.fixture
def sample_rna_sequence():
    """Sample RNA sequence for testing."""
    return "AUCGAUCGAUCG"


@pytest.fixture
def sample_smiles():
    """Sample SMILES string for testing."""
    return "CC(=O)N[C@@H](C)C(=O)O"


@pytest.fixture
def sample_ccd():
    """Sample CCD code for testing."""
    return "ATP"


@pytest.fixture
def protein_chain(sample_protein_sequence):
    """Create a sample protein chain."""
    return ProteinChain(["A"], sample_protein_sequence)


@pytest.fixture
def dna_chain(sample_dna_sequence):
    """Create a sample DNA chain."""
    return DnaChain(["B"], sample_dna_sequence)


@pytest.fixture
def rna_chain(sample_rna_sequence):
    """Create a sample RNA chain."""
    return RnaChain(["C"], sample_rna_sequence)


@pytest.fixture
def ligand_smiles(sample_smiles):
    """Create a sample ligand with SMILES."""
    return Ligand(["L"], smiles=sample_smiles)


@pytest.fixture
def ligand_ccd(sample_ccd):
    """Create a sample ligand with CCD."""
    return Ligand(["M"], ccd=sample_ccd)


@pytest.fixture
def bond():
    """Create a sample bond constraint."""
    return Bond(("A", 1, "CA"), ("B", 2, "CB"))


@pytest.fixture
def contact():
    """Create a sample contact constraint."""
    return Contact(("A", 1), ("B", 2), max_distance=6.0)


@pytest.fixture
def pocket():
    """Create a sample pocket constraint."""
    return Pocket("L", contacts=[("A", 5), ("A", 10)], max_distance=5.0)


@pytest.fixture
def template():
    """Create a sample template."""
    return Template("test.cif", chain_id=["A"], template_id=["X"])


@pytest.fixture
def affinity():
    """Create a sample affinity property."""
    return Affinity("L")


@pytest.fixture
def basic_job(protein_chain, ligand_smiles):
    """Create a basic job with protein and ligand."""
    job = Job("test_job")
    job.sequences = [protein_chain, ligand_smiles]
    return job


@pytest.fixture
def complex_job(
    protein_chain, dna_chain, ligand_smiles, bond, contact, pocket, template, affinity
):
    """Create a complex job with all component types."""
    job = Job("complex_test_job", version=1)
    job.sequences = [protein_chain, dna_chain, ligand_smiles]
    job.constraints = [bond, contact, pocket]
    job.templates = [template]
    job.properties = [affinity]
    return job


@pytest.fixture
def temp_yaml_file():
    """Create a temporary YAML file for testing."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
        yield f.name
    Path(f.name).unlink(missing_ok=True)


@pytest.fixture
def sample_yaml_content():
    """Sample YAML content for testing."""
    return """version: 1
sequences:
  - protein:
      id: A
      sequence: MVTPEGNVSLVDESLLVGVTDED
  - ligand:
      id: L
      smiles: 'CC(=O)N[C@@H](C)C(=O)O'
constraints:
  - pocket:
      binder: L
      contacts: [[A, 5], [A, 10]]
      max_distance: 5.0
properties:
  - affinity:
      binder: L
"""


@pytest.fixture
def invalid_protein_sequence():
    """Invalid protein sequence with illegal characters."""
    return "MVTPXZNB123"


@pytest.fixture
def invalid_dna_sequence():
    """Invalid DNA sequence with illegal characters."""
    return "ATCGXYZ"


@pytest.fixture
def invalid_rna_sequence():
    """Invalid RNA sequence with illegal characters."""
    return "AUCGXYZ"

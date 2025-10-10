"""
This module provides data structures representing various molecular elements, including protein chains, nucleotide
chains (DNA/RNA), ligands/ions, and associated modifications. It is tailored for constructing job configurations
for Boltz-2.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Self, override

from .utils import (
    FlowStyleList,
    SingleQuoted,
    validate_distance,
)

# type definitions
Atom = tuple[str, int, str]  # [CHAIN_ID, RES_IDX, ATOM_NAME]
Token = tuple[str, int | str]  # [CHAIN_ID, RES_IDX/ATOM_NAME]


@dataclass
class Affinity:
    binder: str

    @override
    def __str__(self) -> str:
        return f"++ Affinity for ligand: {self.binder}"

    def to_dict(self) -> dict[str, Any]:
        """Convert the affinity to a dictionary."""
        return {"affinity": {"binder": self.binder}}

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> Self:
        """Create an Affinity instance from a dictionary."""
        return cls(data.get("binder", ""))


@dataclass
class Bond:
    """Represents a bonded atom pair."""

    atom1: Atom
    atom2: Atom

    @override
    def __str__(self) -> str:
        return f"++ Bond between {self.atom1} - {self.atom2}"

    def to_dict(self) -> dict[str, Any]:
        """Convert the bond to a dictionary."""
        return {
            "bond": {
                "atom1": FlowStyleList(self.atom1),
                "atom2": FlowStyleList(self.atom2),
            }
        }

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> Self:
        """Create a Bond instance from a dictionary."""
        atom1 = tuple(data.get("atom1", []))
        atom2 = tuple(data.get("atom2", []))
        return cls(atom1, atom2)


@dataclass
class Pocket:
    """Represents a pocket in the job definition."""

    binder: str
    contacts: list[Token] = field(default_factory=list)
    max_distance: float = 6.0
    force: bool = False

    def __post_init__(self) -> None:
        """Validate distance parameter."""
        validate_distance(self.max_distance)

    @override
    def __str__(self) -> str:
        lines = [
            f"++ Pocket for binder: {self.binder}",
            f"   Max distance: {self.max_distance:.2f} Å",
        ]
        if self.force:
            lines.append("   Force: True")
        if self.contacts:
            lines.append("   Contacts:")
            for chain_id, res_idx_or_atom_name in self.contacts:
                lines.append(f"      - {chain_id} {res_idx_or_atom_name}")
        else:
            lines.append("   No contacts defined.")
        return "\n".join(lines)

    def add_contact_token(self, chain_id: str, res_idx_or_atom_name: int | str) -> None:
        """Add a contact residue/atom."""
        self.contacts.append((chain_id, res_idx_or_atom_name))

    def to_dict(self) -> dict[str, Any]:
        """Convert the pocket to a dictionary."""
        d: dict[str, Any] = {
            "binder": self.binder,
            "contacts": FlowStyleList(self.contacts),
            "max_distance": self.max_distance,
        }
        if self.force:
            d["force"] = self.force
        return {"pocket": d}

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> Self:
        """Create a Pocket instance from a dictionary."""
        binder = data.get("binder", "")
        contacts = data.get("contacts", [])
        max_distance = data.get("max_distance", 6.0)
        force = data.get("force", False)
        return cls(binder, contacts, max_distance, force)


@dataclass
class Contact:
    """Represents a non-bonded contact."""

    token1: Token
    token2: Token
    max_distance: float = 6.0
    force: bool = False

    def __post_init__(self) -> None:
        """Validate distance parameter."""
        validate_distance(self.max_distance)

    @override
    def __str__(self) -> str:
        force_str = ", force: True" if self.force else ""
        return f"++ Contact between {self.token1} - {self.token2}, max distance: {self.max_distance:.2f} Å{force_str}"

    def to_dict(self) -> dict[str, Any]:
        """Convert the contact to a dictionary."""
        d: dict[str, Any] = {
            "token1": FlowStyleList(self.token1),
            "token2": FlowStyleList(self.token2),
            "max_distance": self.max_distance,
        }
        if self.force:
            d["force"] = self.force
        return {"contact": d}

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> Self:
        """Create a Contact instance from a dictionary."""
        token1 = tuple(data.get("token1", []))
        token2 = tuple(data.get("token2", []))
        max_distance = data.get("max_distance", 6.0)
        force = data.get("force", False)
        return cls(token1, token2, max_distance, force)


@dataclass
class SequenceModification:
    """Base class for sequence modifications."""

    ccd: str  # modification type
    position: int  # position of the modification (1-based)

    @override
    def __str__(self) -> str:
        return f"   Modification: {self.ccd} at position {self.position}"

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> Self:
        return cls(
            data.get("ccd", ""),
            data.get("position", 1),
        )


@dataclass
class Template:
    """Represents a template structure for structure prediction.

    Templates can be provided as either mmCIF (.cif) or PDB (.pdb) files.
    The cif and pdb fields are mutually exclusive - exactly one must be provided.

    When force=True, the template structure is enforced during prediction
    using a potential function with the specified threshold distance.
    """

    cif: str | None = None  # mmCIF template file path
    pdb: str | None = None  # PDB template file path
    chain_id: list[str] = field(default_factory=list)
    template_id: list[str] = field(default_factory=list)
    force: bool = False  # Use potential to enforce template
    threshold: float | None = None  # Distance threshold in Angstroms when force=True

    def __post_init__(self) -> None:
        """Validate Template parameters."""
        # Check if both files are provided (mutually exclusive)
        if self.cif and self.pdb:
            raise ValueError(
                "'cif' and 'pdb' are mutually exclusive - provide only one"
            )

        # Exactly one of cif or pdb must be provided (and not empty)
        if not self.cif and not self.pdb:
            raise ValueError("Either 'cif' or 'pdb' file path must be provided")

        # If force is True, threshold must be specified
        if self.force and self.threshold is None:
            raise ValueError("'threshold' must be specified when 'force=True'")

    @override
    def __str__(self) -> str:
        # Determine file type and path
        if self.cif:
            file_info = f"CIF: {self.cif!r}"
        elif self.pdb:
            file_info = f"PDB: {self.pdb!r}"
        else:
            file_info = "<no file specified>"

        lines = [f"++ Template: {file_info}"]

        if self.template_id:
            lines.append(f"   Template ID(s): {', '.join(self.template_id)}")
        if self.chain_id:
            lines.append(f"   For chain ID(s): {', '.join(self.chain_id)}")
        if self.force:
            lines.append(f"   Force: True (threshold: {self.threshold}Å)")

        return "\n".join(lines)

    def to_dict(self) -> dict[str, Any]:
        """Convert the template to a dictionary."""
        d: dict[str, Any] = {}

        # Add the file path (either cif or pdb)
        if self.cif:
            d["cif"] = self.cif
        elif self.pdb:
            d["pdb"] = self.pdb

        if self.chain_id:
            d["chain_id"] = (
                FlowStyleList(self.chain_id)
                if len(self.chain_id) > 1
                else self.chain_id[0]
            )

        if self.template_id:
            d["template_id"] = (
                FlowStyleList(self.template_id)
                if len(self.template_id) > 1
                else self.template_id[0]
            )

        # Add force and threshold if applicable
        if self.force:
            d["force"] = self.force
        if self.threshold is not None:
            d["threshold"] = self.threshold

        return d

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> Self:
        """Create a Template instance from a dictionary."""
        cif = data.get("cif")
        pdb = data.get("pdb")

        if chain_id := data.get("chain_id"):
            if isinstance(chain_id, str):
                chain_id = [chain_id]
        else:
            chain_id = []

        if template_id := data.get("template_id"):
            if isinstance(template_id, str):
                template_id = [template_id]
        else:
            template_id = []

        force = data.get("force", False)
        threshold = data.get("threshold")

        return cls(
            cif=cif,
            pdb=pdb,
            chain_id=chain_id,
            template_id=template_id,
            force=force,
            threshold=threshold,
        )


@dataclass
class Chain:
    """Base class for protein and nucleotide chains."""

    ids: list[str]
    sequence: str
    modifications: list[SequenceModification] = field(default_factory=list)
    cyclic: bool = False

    def __post_init__(self) -> None:
        """Check for empty sequence."""
        if not self.sequence:
            raise ValueError("Sequence cannot be empty.")

    @override
    def __str__(self) -> str:
        if len(self.sequence) > 25:
            seq = f"{self.sequence[:10]}.....{self.sequence[-10:]}"
        else:
            seq = self.sequence
        lines = [
            f"   {'ID' if len(self.ids) == 1 else 'IDs'}: {', '.join(self.ids)}",
            f"   Sequence: {seq}, {len(self.sequence)} residues",
        ]
        if self.cyclic:
            lines.append("   Cyclic: True")
        lines.extend(str(m) for m in self.modifications)
        return "\n".join(lines)

    def add_modification(self, ccd: str, position: int) -> Self:
        """Add a modification to the chain."""
        mod = SequenceModification(ccd, position)
        self.modifications.append(mod)
        return self

    def set_cyclic(self, cyclic: bool) -> None:
        """Set the cyclic property of the chain."""
        if not isinstance(cyclic, bool):
            raise TypeError("Cyclic must be a boolean value.")
        self.cyclic = cyclic

    def to_dict(self) -> dict[str, Any]:
        """Convert the chain to a dictionary."""
        d: dict[str, Any] = {
            "id": FlowStyleList(self.ids) if len(self.ids) > 1 else self.ids[0],
            "sequence": self.sequence,
        }

        if self.modifications:
            d["modifications"] = [
                {"position": mod.position, "ccd": mod.ccd} for mod in self.modifications
            ]

        if self.cyclic:
            d["cyclic"] = self.cyclic

        return d

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> Self:
        ids = data.get("id", [])
        if isinstance(ids, str):
            ids = [ids]
        chain = cls(
            ids=ids,
            sequence=data.get("sequence", ""),
            cyclic=data.get("cyclic", False),
        )

        for mod in data.get("modifications", []):
            chain.modifications.append(SequenceModification.from_dict(mod))

        return chain


@dataclass
class ProteinChain(Chain):
    """Represents a protein chain in the job definition."""

    msa: str | None = None

    def __post_init__(self) -> None:
        """Check if the protein sequence contains only legal one-letter codes."""
        super().__post_init__()
        one_letter_codes = set("ACDEFGHIKLMNPQRSTVWY")
        diff = set(self.sequence.upper()).difference(one_letter_codes)
        if diff:
            raise ValueError(
                f"Protein sequence contains invalid one-letter codes: {', '.join(repr(b) for b in sorted(diff))}."
            )

    @override
    def __str__(self) -> str:
        lines = [
            "++ Protein chain:",
            super().__str__(),
        ]

        if self.msa:
            lines.append(f"   MSA: {self.msa!r}")

        return "\n".join(lines)

    def set_msa_path(self, path: str) -> Self:
        """Set the path to the MSA file.

        Args:
            path: Path to the MSA file or special mode selector:
                - `.a3m` file: Standard MSA format for single protein chains
                - `.csv` file: Multi-chain MSA format with two columns:
                  - `sequence`: protein sequence
                  - `key`: unique identifier for matching rows across chains
                  (sequences with the same key are mutually aligned)
                - "empty": Force single-sequence mode (not recommended, reduces accuracy)

        Returns:
            Self for method chaining

        Note:
            For multi-chain complexes, use CSV format to ensure proper MSA pairing
            across protein chains via the `key` column matching system.
        """
        self.msa = path
        return self

    @override
    def to_dict(self) -> dict[str, Any]:
        """Convert the protein chain to a dictionary."""
        d = super().to_dict()
        if self.msa is not None:
            d["msa"] = self.msa

        return {"protein": d}

    @override
    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> Self:
        """Create a ProteinChain instance from a dictionary."""
        chain = super().from_dict(data)
        return cls(
            ids=chain.ids,
            sequence=chain.sequence,
            modifications=chain.modifications,
            cyclic=chain.cyclic,
            msa=data.get("msa"),
        )


class DnaChain(Chain):
    """Represents a DNA chain in the job definition."""

    def __post_init__(self) -> None:
        """Check if the DNA sequence contains only 'A', 'C', 'G', or 'T'."""
        super().__post_init__()
        diff = set(self.sequence.upper()).difference({"A", "C", "G", "T"})
        if diff:
            raise ValueError(
                f"DNA sequence can only contain 'A', 'C', 'G', or 'T'. Found {', '.join(repr(b) for b in sorted(diff))}."
            )

    @override
    def __str__(self) -> str:
        lines = [
            "++ DNA chain:",
            super().__str__(),
        ]
        return "\n".join(lines)

    @override
    def to_dict(self) -> dict[str, Any]:
        """Convert the DNA chain to a dictionary."""
        d = super().to_dict()
        return {"dna": d}

    @override
    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> Self:
        """Create a ProteinChain instance from a dictionary."""
        chain = super().from_dict(data)
        return cls(
            ids=chain.ids,
            sequence=chain.sequence,
            modifications=chain.modifications,
            cyclic=chain.cyclic,
        )


class RnaChain(Chain):
    """Represents an RNA chain in the job definition."""

    def __post_init__(self) -> None:
        """Check if the RNA sequence contains only 'A', 'C', 'G', or 'U'."""
        super().__post_init__()
        diff = set(self.sequence.upper()).difference({"A", "C", "G", "U"})
        if diff:
            raise ValueError(
                f"RNA sequence can only contain 'A', 'C', 'G', or 'U'. Found {', '.join(repr(b) for b in sorted(diff))}."
            )

    @override
    def __str__(self) -> str:
        lines = [
            "++ RNA chain:",
            super().__str__(),
        ]
        return "\n".join(lines)

    @override
    def to_dict(self) -> dict[str, Any]:
        """Convert the RNA chain to a dictionary."""
        d = super().to_dict()
        return {"rna": d}

    @override
    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> Self:
        """Create a ProteinChain instance from a dictionary."""
        chain = super().from_dict(data)
        return cls(
            ids=chain.ids,
            sequence=chain.sequence,
            modifications=chain.modifications,
            cyclic=chain.cyclic,
        )


@dataclass
class Ligand:
    """Represents a ligand or an ion in the job definition."""

    ids: list[str]
    ccd: None | str = None
    smiles: str = ""

    @override
    def __str__(self) -> str:
        lines = [
            "++ Ligand:",
            f"   {'ID' if len(self.ids) == 1 else 'IDs'}: {', '.join(self.ids)}",
        ]
        if self.smiles:
            lines.append(f"   SMILES: {self.smiles}")
        elif self.ccd:
            lines.append(f"   CCD: {self.ccd}")
        return "\n".join(lines)

    def to_dict(self) -> dict[str, Any]:
        """Convert the ligand to a dictionary."""
        d: dict[str, Any] = {
            "id": FlowStyleList(self.ids) if len(self.ids) > 1 else self.ids[0],
        }

        # if a SMILES string is provided, use it; otherwise, use CCD codes
        if self.smiles:
            d["smiles"] = SingleQuoted(self.smiles)
        elif self.ccd:
            d["ccd"] = self.ccd
        else:
            # This should not happen, if the ligand object is created via the Job class
            raise RuntimeError(
                "Either a SMILES or a CCD code must be provided for the ligand."
            )
        return {"ligand": d}

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> Self:
        ids = data.get("id")
        if isinstance(ids, str):
            ids = [ids]
        return cls(ids, data.get("ccd"), data.get("smiles") or "")  # type: ignore

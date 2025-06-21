"""
This module provides data structures representing various molecular elements, including protein chains, nucleotide
chains (DNA/RNA), ligands/ions, and associated modifications. It is tailored for constructing job configurations
for Boltz-2.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Self

from .utils import (
    FlowStyleList,
    SingleQuoted,
)

# type definitions
Atom = tuple[str, int, str]  # [CHAIN_ID, RES_IDX, ATOM_NAME]
Token = tuple[str, int | str]  # [CHAIN_ID, RES_IDX/ATOM_NAME]


@dataclass
class Affinity:
    binder: str

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
    max_distance: float = 5.0

    def __str__(self) -> str:
        lines = [
            f"++ Pocket for binder: {self.binder}",
            f"   Max distance: {self.max_distance:.2f} Å",
        ]
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
        return {
            "pocket": {
                "binder": self.binder,
                "contacts": FlowStyleList(self.contacts),
                "max_distance": self.max_distance,
            }
        }

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> Self:
        """Create a Pocket instance from a dictionary."""
        binder = data.get("binder", "")
        contacts = data.get("contacts", [])
        max_distance = data.get("max_distance", 5.0)
        return cls(binder, contacts, max_distance)


@dataclass
class Contact:
    """Represents a non-bonded contact."""

    token1: Token
    token2: Token
    max_distance: float = 5.0

    def __str__(self) -> str:
        return f"++ Contact between {self.token1} - {self.token2}, max distance: {self.max_distance:.2f} Å"

    def to_dict(self) -> dict[str, Any]:
        """Convert the contact to a dictionary."""
        return {
            "contact": {
                "token1": FlowStyleList(self.token1),
                "token2": FlowStyleList(self.token2),
                "max_distance": self.max_distance,
            }
        }

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> Self:
        """Create a Contact instance from a dictionary."""
        token1 = tuple(data.get("token1", []))
        token2 = tuple(data.get("token2", []))
        max_distance = data.get("max_distance", 5.0)
        return cls(token1, token2, max_distance)


@dataclass
class SequenceModification:
    """Base class for sequence modifications."""

    mod_type: str  # modification type
    position: int  # position of the modification (1-based)

    def __str__(self) -> str:
        return f"   Modification: {self.mod_type} at position {self.position}"

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> Self:
        return cls(
            data.get("mad_type"),  # type: ignore
            data.get("position"),  # type: ignore
        )


@dataclass
class Template:
    """Represents a template for a protein chain."""

    cif: str  # mmCIF template string
    chain_id: list[str] = field(default_factory=list)
    template_id: list[str] = field(default_factory=list)

    def __str__(self) -> str:
        lines = [
            f"    {self.cif!r}",
        ]
        if self.template_id:
            lines.append(f"template id(s): {', '.join(self.template_id)}")
        if self.chain_id:
            lines.append(f"to be used for chain id(s): {', '.join(self.chain_id)}")
        return ", ".join(lines)

    def to_dict(self) -> dict[str, Any]:
        """Convert the template to a dictionary."""
        d: dict[str, Any] = {
            "cif": self.cif,
        }
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

        return d

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> Self:
        cif = data.get("cif")
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
        return cls(cif, chain_id, template_id)  # type: ignore


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

    def add_modification(self, mod_type: str, position: int) -> Self:
        """Add a modification to the chain."""
        mod = SequenceModification(mod_type, position)
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
                {"position": mod.position, "ccd": mod.mod_type}
                for mod in self.modifications
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

    def __str__(self) -> str:
        lines = [
            "++ Protein chain:",
            super().__str__(),
        ]

        if self.msa:
            lines.append(f"   MSA: {self.msa!r}")

        return "\n".join(lines)

    def set_msa_path(self, path: str) -> Self:
        """Set the path to the unpaired MSA file. ONLY for AF3 input file version >= 2."""
        self.msa = path
        return self

    def to_dict(self) -> dict[str, Any]:
        """Convert the protein chain to a dictionary."""
        d = super().to_dict()
        if self.msa is not None:
            d["msa"] = self.msa

        return {"protein": d}

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

    def __str__(self) -> str:
        lines = [
            "++ DNA chain:",
            super().__str__(),
        ]
        return "\n".join(lines)

    def to_dict(self) -> dict[str, Any]:
        """Convert the DNA chain to a dictionary."""
        d = super().to_dict()
        return {"dna": d}

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

    def __str__(self) -> str:
        lines = [
            "++ RNA chain:",
            super().__str__(),
        ]
        return "\n".join(lines)

    def to_dict(self) -> dict[str, Any]:
        """Convert the RNA chain to a dictionary."""
        d = super().to_dict()
        return {"rna": d}

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

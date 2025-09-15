"""
Module for defining molecular components, modifications, and job configurations for Boltz-2.

The `Job` class is the main container for combining chains, ligands/ions, and modifications that can be converted
to YAML as input for Boltz-2.

It follows the YAML schema defined in the Boltz-2 documentation (as of June 2024):
https://github.com/jwohlwend/boltz/blob/main/docs/prediction.md
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field
from itertools import islice
from typing import Any, Self, Type

import yaml

from .components import (
    Affinity,
    Bond,
    Contact,
    DnaChain,
    Ligand,
    Pocket,
    ProteinChain,
    RnaChain,
    Template,
)
from .utils import chain_id, IndentedDumper

# type definitions
Sequence = ProteinChain | DnaChain | RnaChain | Ligand
Constraint = Bond | Contact | Pocket
Property = Affinity


@dataclass
class Job:
    """Represents an Boltz-2 job with methods to add entities."""

    name: str = "boltz_job"  # default job name
    version: int = 1  # default version
    sequences: list[Sequence] = field(default_factory=list)
    constraints: list[Constraint] = field(default_factory=list)
    properties: list[Affinity] = field(default_factory=list)
    templates: list[Template] = field(default_factory=list)

    def __post_init__(self) -> None:
        """Post-initialization method."""
        # Initialize the chain ID generator
        self._chain_ids = chain_id()

    def __str__(self) -> str:
        """Return a string representation of the job."""
        lines = [
            f"Job name: {self.name}",
            f"Version: {self.version}",
        ]
        num_seq = len(self.sequences)
        if num_seq == 0:
            lines.append("No sequences.")
        else:
            lines.append(f"{num_seq} Sequence(s):")
            for seq in self.sequences:
                lines.append(str(seq))

        if self.constraints:
            lines.append("Constraints:")
            for constraint in self.constraints:
                lines.append(str(constraint))

        if self.templates:
            lines.append("Templates:")
            for template in self.templates:
                lines.append(str(template))

        if self.properties:
            lines.append("Properties:")
            for property in self.properties:
                lines.append(str(property))

        return "\n".join(lines)

    def _get_current_ids(self) -> list[str]:
        """Get the current chain IDs."""
        seq_ids: list[str] = []
        for seq in self.sequences:
            seq_ids.extend(seq.ids)
        return seq_ids

    def _check_ids(self) -> None:
        """Check for duplicate IDs."""
        seq_ids = self._get_current_ids()
        if len(seq_ids) != len(set(seq_ids)):
            raise ValueError(f"Duplicate chain IDs found: {seq_ids}")

    def _get_ids(self, ids: None | str | list[str], count: int) -> list[str]:
        if ids is None:
            ids = list(islice(self._chain_ids, count))
        else:
            if isinstance(ids, str):
                ids = [ids]
            elif not isinstance(ids, list):
                raise TypeError("IDs must be a string or a list of strings.")
        if not ids:
            raise ValueError("Number of chains or ligands must be greater than zero.")
        return ids

    def add_protein_chain(
        self,
        sequence: str,
        count: int = 1,
        ids: None | str | list[str] = None,
        cyclic: bool = False,
        msa: str | None = None,
    ) -> ProteinChain:
        """Add a protein chain to the job."""
        chn = ProteinChain(self._get_ids(ids, count), sequence, cyclic=cyclic, msa=msa)
        self.sequences.append(chn)
        return chn

    def add_dna_chain(
        self,
        sequence: str,
        count: int = 1,
        ids: None | str | list[str] = None,
        cyclic: bool = False,
    ) -> DnaChain:
        """Add a DNA chain to the job."""
        chn = DnaChain(self._get_ids(ids, count), sequence, cyclic=cyclic)
        self.sequences.append(chn)
        return chn

    def add_rna_chain(
        self,
        sequence: str,
        count: int = 1,
        ids: None | str | list[str] = None,
        cyclic: bool = False,
    ) -> RnaChain:
        """Add an RNA chain to the job."""
        chn = RnaChain(self._get_ids(ids, count), sequence, cyclic=cyclic)
        self.sequences.append(chn)
        return chn

    def add_ligand(
        self,
        ccd: str | None = None,
        smiles: str = "",
        count: int = 1,
        ids: None | str | list[str] = None,
    ) -> Ligand:
        """Add a ligand to the job."""
        if not ccd and not smiles:
            raise ValueError("Either CCD codes or SMILES string must be provided.")
        if ccd and smiles:
            warnings.warn(
                "`ccd` and `smiles` are given - they are mutually exclusive - will be using smiles"
            )
        _ids = self._get_ids(ids, count)
        if smiles:
            ligand = Ligand(_ids, smiles=smiles)
        else:
            ligand = Ligand(_ids, ccd=ccd)
        self.sequences.append(ligand)
        return ligand

    def add_bond(
        self, id_1: str, resi_1: int, name_1: str, id_2: str, resi_2: int, name_2: str
    ) -> None:
        """Add a bond constraint."""
        seg_ids = self._get_current_ids()
        if id_1 not in seg_ids or id_2 not in seg_ids:
            raise ValueError(
                f"Both chain IDs must be defined: {id_1}, {id_2}. "
                f"Current IDs: {seg_ids}"
            )
        self.constraints.append(Bond((id_1, resi_1, name_1), (id_2, resi_2, name_2)))

    def add_disulfide_bond(
        self, id_1: str, resi_1: int, id_2: str, resi_2: int
    ) -> None:
        """Add a disulfide bond constraint."""
        seg_ids = self._get_current_ids()
        if id_1 not in seg_ids or id_2 not in seg_ids:
            raise ValueError(
                f"Both chain IDs must be defined: {id_1}, {id_2}. "
                f"Current IDs: {seg_ids}"
            )
        # check that both residues are cysteines
        for seq in self.sequences:
            if id_1 in seq.ids and isinstance(seq, ProteinChain):
                if seq.sequence[resi_1 - 1] != "C":
                    raise ValueError(
                        f"Residue {resi_1} in chain {id_1} is not a cysteine."
                    )
            if id_2 in seq.ids and isinstance(seq, ProteinChain):
                if seq.sequence[resi_2 - 1] != "C":
                    raise ValueError(
                        f"Residue {resi_2} in chain {id_2} is not a cysteine."
                    )
        self.constraints.append(Bond((id_1, resi_1, "SG"), (id_2, resi_2, "SG")))

    def add_pocket(
        self,
        binder: str,
        contact_tokens: None | list[tuple[str, int | str]] = None,
        max_distance: float = 5.0,
    ) -> Pocket:
        """Add a pocket constraint for a ligand with the given binder ID."""
        if contact_tokens is None:
            contact_tokens = []
        for seq in self.sequences:
            if binder in seq.ids and isinstance(seq, Ligand):
                pocket = Pocket(
                    binder, contacts=contact_tokens, max_distance=max_distance
                )
                self.constraints.append(pocket)
                return pocket
        else:
            raise ValueError(
                f"Pocket can only be defined for ligands. No ligand with id {binder} found."
            )

    def add_contact(
        self,
        id_1: str,
        resi_1: int | str,
        id_2: str,
        resi_2: int | str,
        max_distance: float = 5.0,
    ) -> None:
        """Add a bond constraint."""
        seg_ids = self._get_current_ids()
        if id_1 not in seg_ids or id_2 not in seg_ids:
            raise ValueError(
                f"Both chain IDs must be defined: {id_1}, {id_2}. "
                f"Current IDs: {seg_ids}"
            )
        self.constraints.append(
            Contact((id_1, resi_1), (id_2, resi_2), max_distance=max_distance)
        )

    def add_template(
        self,
        cif: str,
        chain_id: str | list[str] | None = None,
        template_id: str | list[str] | None = None,
    ) -> None:
        """Add a template structure."""
        if not cif:
            raise ValueError("CIF file must be provided.")
        if chain_id is None:
            chain_id = []
        elif isinstance(chain_id, str):
            chain_id = [chain_id]
        if template_id is None:
            template_id = []
        elif isinstance(template_id, str):
            template_id = [template_id]

        self.templates.append(Template(cif, chain_id=chain_id, template_id=template_id))

    def request_affinity(self, binder: str) -> None:
        """Request an affinity estimation for a ligand with the given binder ID."""
        for seq in self.sequences:
            if binder in seq.ids and isinstance(seq, Ligand):
                self.properties.append(Affinity(binder))
                break
        else:
            raise ValueError(
                f"Affinity can only be estimated for ligands. No ligand with id {binder} found."
            )

    def to_dict(self) -> dict[str, Any]:
        """Convert the Job to a dictionary suitable for JSON serialization."""
        d: dict[str, Any] = {
            "version": self.version,
        }

        # add sequences / ligands / ions
        if self.sequences:
            self._check_ids()
            d["sequences"] = [seq.to_dict() for seq in self.sequences]
        else:
            raise ValueError("Empty list of sequences.")

        if self.constraints:
            d["constraints"] = [constraint.to_dict() for constraint in self.constraints]

        if self.templates:
            d["templates"] = [template.to_dict() for template in self.templates]

        if self.properties:
            d["properties"] = [prop.to_dict() for prop in self.properties]

        return d

    def write_yaml(self, filename: str) -> None:
        """Write the job to a YAML file as input for Boltz-2."""
        with open(filename, "w") as f:
            yaml.dump(
                self.to_dict(), f, Dumper=IndentedDumper, indent=2, sort_keys=False
            )

    @classmethod
    def from_yaml(cls, filename: str, jobname: str = "boltz_job") -> Self:
        """Read a job from a YAML file."""
        sequence_match: dict[str, Type[Sequence]] = {
            "protein": ProteinChain,
            "dna": DnaChain,
            "rna": RnaChain,
            "ligand": Ligand,
        }

        constraint_match: dict[str, Type[Constraint]] = {
            "bond": Bond,
            "contact": Contact,
            "pocket": Pocket,
        }

        property_match: dict[str, Type[Property]] = {
            "affinity": Affinity,
        }

        with open(filename) as f:
            data = yaml.safe_load(f)
        job = cls(jobname, version=data.get("version", 1))

        if sequences := data.get("sequences"):
            for sequence in sequences:
                ((seq_type, seq_data),) = sequence.items()
                job.sequences.append(sequence_match[seq_type].from_dict(seq_data))

        if constraints := data.get("constraints"):
            for constraint in constraints:
                ((const_type, const_data),) = constraint.items()
                job.constraints.append(
                    constraint_match[const_type].from_dict(const_data)
                )

        if templates := data.get("templates"):
            for template in templates:
                job.templates.append(Template.from_dict(template))

        if properties := data.get("properties"):
            for property in properties:
                ((prop_type, prop_data),) = property.items()
                job.properties.append(property_match[prop_type].from_dict(prop_data))

        return job

"""
Module for defining molecular components, modifications, and job configurations for Boltz-2.

The `Job` class is the main container for combining chains, ligands/ions, and modifications that can be converted
to YAML as input for Boltz-2.

It follows the YAML schema defined in the Boltz-2 documentation (as of June 2024):
https://github.com/jwohlwend/boltz/blob/main/docs/prediction.md
"""

from __future__ import annotations

import base64
import warnings
from dataclasses import dataclass, field
from itertools import islice
from pathlib import Path
from typing import Any, Self, override
from urllib.parse import urlsplit

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

ApiDict = dict[str, Any]


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

    @override
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
                f"Both chain IDs must be defined: {id_1}, {id_2}. Current IDs: {seg_ids}"
            )
        self.constraints.append(Bond((id_1, resi_1, name_1), (id_2, resi_2, name_2)))

    def add_disulfide_bond(
        self, id_1: str, resi_1: int, id_2: str, resi_2: int
    ) -> None:
        """Add a disulfide bond constraint."""
        # Check that both residues are cysteines
        chain_1_found = False
        chain_2_found = False
        for seq in self.sequences:
            if isinstance(seq, ProteinChain):
                if id_1 in seq.ids:
                    chain_1_found = True
                    if seq.sequence[resi_1 - 1] != "C":
                        raise ValueError(
                            f"Residue {resi_1} in chain {id_1} is not a cysteine."
                        )
                if id_2 in seq.ids:
                    chain_2_found = True
                    if seq.sequence[resi_2 - 1] != "C":
                        raise ValueError(
                            f"Residue {resi_2} in chain {id_2} is not a cysteine."
                        )
        if not chain_1_found or not chain_2_found:
            raise ValueError("One or both chain IDs not found as ProteinChain.")
        if id_1 == id_2 and resi_1 == resi_2:
            raise ValueError(
                f"Cannot form a disulfide bond with the same residue. {id_1}-{resi_1}, {id_2}-{resi_2}"
            )

        self.constraints.append(Bond((id_1, resi_1, "SG"), (id_2, resi_2, "SG")))

    def add_pocket(
        self,
        binder: str,
        contact_tokens: None | list[tuple[str, int | str]] = None,
        max_distance: float = 6.0,
        force: bool = False,
    ) -> Pocket:
        """Add a pocket constraint for any sequence with the given binder ID.

        The binder can be a molecule, protein, DNA or RNA as per Boltz-2 schema.
        """
        if contact_tokens is None:
            contact_tokens = []

        # Check if binder ID exists in any sequence
        seg_ids = self._get_current_ids()
        if binder not in seg_ids:
            raise ValueError(
                f"No sequence with id '{binder}' found. Available IDs: {seg_ids}"
            )

        pocket = Pocket(
            binder, contacts=contact_tokens, max_distance=max_distance, force=force
        )
        self.constraints.append(pocket)
        return pocket

    def add_contact(
        self,
        id_1: str,
        resi_1: int | str,
        id_2: str,
        resi_2: int | str,
        max_distance: float = 6.0,
        force: bool = False,
    ) -> Contact:
        """Add a contact constraint between two sequence positions."""
        seg_ids = self._get_current_ids()
        if id_1 not in seg_ids or id_2 not in seg_ids:
            raise ValueError(
                f"Both chain IDs must be defined: {id_1}, {id_2}. Current IDs: {seg_ids}"
            )
        contact = Contact(
            (id_1, resi_1), (id_2, resi_2), max_distance=max_distance, force=force
        )
        self.constraints.append(contact)
        return contact

    def add_template(
        self,
        cif: str | None = None,
        pdb: str | None = None,
        chain_id: str | list[str] | None = None,
        template_id: str | list[str] | None = None,
        force: bool = False,
        threshold: float | None = None,
    ) -> Template:
        """Add a template structure from a CIF or PDB file.

        Args:
            cif: Path to mmCIF template file (mutually exclusive with pdb)
            pdb: Path to PDB template file (mutually exclusive with cif)
            chain_id: Chain ID(s) in the job to apply this template to
            template_id: Template chain ID(s) from the template file
            force: Whether to use potential to enforce template structure
            threshold: Distance threshold in Angstroms when force=True

        Returns:
            The created Template instance

        Raises:
            ValueError: If neither cif nor pdb is provided, or if both are provided,
                       or if force=True but threshold is None
        """
        # Backward compatibility: if first argument is a string, treat as cif
        if isinstance(cif, str) and pdb is None:
            pass  # Normal case
        elif cif is None and isinstance(pdb, str):
            pass  # PDB case
        elif cif is None and pdb is None:
            raise ValueError("Either 'cif' or 'pdb' file path must be provided")
        else:
            raise ValueError(
                "'cif' and 'pdb' are mutually exclusive - provide only one"
            )

        # Normalize chain_id and template_id to lists
        if chain_id is None:
            chain_id = []
        elif isinstance(chain_id, str):
            chain_id = [chain_id]

        if template_id is None:
            template_id = []
        elif isinstance(template_id, str):
            template_id = [template_id]

        template = Template(
            cif=cif,
            pdb=pdb,
            chain_id=chain_id,
            template_id=template_id,
            force=force,
            threshold=threshold,
        )
        self.templates.append(template)
        return template

    def add_pdb_template(
        self,
        pdb: str,
        chain_id: str | list[str] | None = None,
        template_id: str | list[str] | None = None,
        force: bool = False,
        threshold: float | None = None,
    ) -> Template:
        """Add a PDB template structure.

        Convenience method for adding PDB templates without specifying cif=None.

        Args:
            pdb: Path to PDB template file
            chain_id: Chain ID(s) in the job to apply this template to
            template_id: Template chain ID(s) from the PDB file
            force: Whether to use potential to enforce template structure
            threshold: Distance threshold in Angstroms when force=True

        Returns:
            The created Template instance
        """
        return self.add_template(
            pdb=pdb,
            chain_id=chain_id,
            template_id=template_id,
            force=force,
            threshold=threshold,
        )

    def add_cif_template(
        self,
        cif: str,
        chain_id: str | list[str] | None = None,
        template_id: str | list[str] | None = None,
        force: bool = False,
        threshold: float | None = None,
    ) -> Template:
        """Add a CIF template structure.

        Convenience method for adding CIF templates without specifying pdb=None.

        Args:
            cif: Path to mmCIF template file
            chain_id: Chain ID(s) in the job to apply this template to
            template_id: Template chain ID(s) from the CIF file
            force: Whether to use potential to enforce template structure
            threshold: Distance threshold in Angstroms when force=True

        Returns:
            The created Template instance
        """
        return self.add_template(
            cif=cif,
            chain_id=chain_id,
            template_id=template_id,
            force=force,
            threshold=threshold,
        )

    def request_affinity(self, binder: str) -> None:
        """Request an affinity estimation for a ligand with the given binder ID.

        Args:
            binder: Chain ID of the ligand to compute affinity for

        Raises:
            ValueError: If binder is not a ligand, if affinity already requested,
                       or if no ligand with the specified ID is found

        Warnings:
            Issues warnings if non-protein targets are present, as this may
            produce unreliable affinity results according to Boltz-2 schema
        """
        # Check if affinity computation has already been requested
        existing_affinities = [
            prop for prop in self.properties if isinstance(prop, Affinity)
        ]
        if existing_affinities:
            raise ValueError(
                f"Only one affinity computation allowed per job. "
                f"Affinity already requested for: {existing_affinities[0].binder}"
            )

        # Find the ligand with the specified binder ID
        binder_ligand = None
        for seq in self.sequences:
            if binder in seq.ids and isinstance(seq, Ligand):
                binder_ligand = seq
                break

        if binder_ligand is None:
            raise ValueError(
                f"Affinity can only be estimated for ligands. No ligand with id '{binder}' found."
            )

        # Check for protein targets and issue warnings for non-protein targets
        non_ligand_sequences = [
            seq for seq in self.sequences if not isinstance(seq, Ligand)
        ]
        protein_targets = [
            seq for seq in non_ligand_sequences if isinstance(seq, ProteinChain)
        ]
        non_protein_targets = [
            seq for seq in non_ligand_sequences if not isinstance(seq, ProteinChain)
        ]

        if not protein_targets:
            warnings.warn(
                "No protein chains found in job. Affinity computation without protein "
                "targets may produce unreliable results.",
                UserWarning,
                stacklevel=2,
            )

        if non_protein_targets:
            target_types = [type(seq).__name__ for seq in non_protein_targets]
            warnings.warn(
                f"Non-protein targets detected: {', '.join(target_types)}. "
                "Affinity computation with DNA/RNA targets may produce unreliable results "
                "according to Boltz-2 schema.",
                UserWarning,
                stacklevel=2,
            )

        # Add the affinity property
        self.properties.append(Affinity(binder))

    def _get_sequence_by_id(self, chain_id: str) -> Sequence:
        """Return the sequence containing the given chain ID."""
        for seq in self.sequences:
            if chain_id in seq.ids:
                return seq
        raise ValueError(f"No sequence with id '{chain_id}' found.")

    @staticmethod
    def _is_url(value: str) -> bool:
        """Return whether a string is an HTTPS URL supported by the Boltz API."""
        return value.startswith("https://")

    @staticmethod
    def _suffix_from_source(source: str) -> str:
        """Return a lowercase suffix from a local path or URL path."""
        if Job._is_url(source):
            return Path(urlsplit(source).path).suffix.lower()
        return Path(source).expanduser().suffix.lower()

    @staticmethod
    def _base64_source(path: str, media_type: str) -> ApiDict:
        """Read a local file and return a Boltz API base64 source object."""
        file_path = Path(path).expanduser()
        if not file_path.is_file():
            raise ValueError(f"Local file does not exist: {path}")
        return {
            "type": "base64",
            "data": base64.b64encode(file_path.read_bytes()).decode("ascii"),
            "media_type": media_type,
        }

    @staticmethod
    def _file_source(path_or_url: str, media_types: dict[str, str]) -> ApiDict:
        """Return a URL or base64 source object for an API file field."""
        suffix = Job._suffix_from_source(path_or_url)
        if suffix not in media_types:
            supported = ", ".join(sorted(media_types))
            raise ValueError(
                f"Unsupported file extension for Boltz API source: {path_or_url!r}. "
                f"Supported extensions: {supported}"
            )
        if Job._is_url(path_or_url):
            return {"type": "url", "url": path_or_url}
        return Job._base64_source(path_or_url, media_types[suffix])

    @staticmethod
    def _api_modifications(seq: ProteinChain | DnaChain | RnaChain) -> list[ApiDict]:
        """Convert 1-based local modifications to 0-based API modifications."""
        modifications = []
        for mod in seq.modifications:
            if mod.position < 1:
                raise ValueError(
                    f"Modification positions are 1-based and must be greater than zero: {mod.position}"
                )
            modifications.append(
                {
                    "residue_index": mod.position - 1,
                    "type": "ccd",
                    "value": mod.ccd,
                }
            )
        return modifications

    @staticmethod
    def _api_msa(msa: str) -> ApiDict:
        """Convert local protein MSA settings to the Boltz API shape."""
        if msa == "empty":
            return {"type": "empty"}

        suffix = Job._suffix_from_source(msa)
        formats = {".a3m": "a3m", ".csv": "csv"}
        if suffix not in formats:
            raise ValueError(
                f"Unsupported MSA extension for Boltz API source: {msa!r}. "
                "Supported extensions: .a3m, .csv"
            )

        return {
            "type": "custom",
            "format": formats[suffix],
            "source": Job._file_source(msa, {".a3m": "text/x-a3m", ".csv": "text/csv"}),
        }

    def _api_entity(self, seq: Sequence) -> ApiDict:
        """Convert a package sequence object to a Boltz API entity."""
        if isinstance(seq, ProteinChain):
            entity: ApiDict = {
                "type": "protein",
                "value": seq.sequence,
                "chain_ids": seq.ids,
            }
            if seq.cyclic:
                entity["cyclic"] = seq.cyclic
            if modifications := self._api_modifications(seq):
                entity["modifications"] = modifications
            if seq.msa is not None:
                entity["msa"] = self._api_msa(seq.msa)
            return entity

        if isinstance(seq, DnaChain):
            entity = {
                "type": "dna",
                "value": seq.sequence,
                "chain_ids": seq.ids,
            }
            if seq.cyclic:
                entity["cyclic"] = seq.cyclic
            if modifications := self._api_modifications(seq):
                entity["modifications"] = modifications
            return entity

        if isinstance(seq, RnaChain):
            entity = {
                "type": "rna",
                "value": seq.sequence,
                "chain_ids": seq.ids,
            }
            if seq.cyclic:
                entity["cyclic"] = seq.cyclic
            if modifications := self._api_modifications(seq):
                entity["modifications"] = modifications
            return entity

        if isinstance(seq, Ligand):
            if seq.smiles:
                return {
                    "type": "ligand_smiles",
                    "value": seq.smiles,
                    "chain_ids": seq.ids,
                }
            if seq.ccd:
                return {
                    "type": "ligand_ccd",
                    "value": seq.ccd,
                    "chain_ids": seq.ids,
                }
            raise ValueError("Ligand must define either a SMILES string or CCD code.")

        raise TypeError(f"Unsupported sequence type: {type(seq).__name__}")

    def _api_polymer_atom(self, chain_id: str, residue_index: int, atom_name: str) -> ApiDict:
        """Convert a package atom tuple to a Boltz API atom reference."""
        if residue_index < 1:
            raise ValueError(
                f"Residue indices are 1-based and must be greater than zero: {residue_index}"
            )
        return {
            "type": "polymer_atom",
            "chain_id": chain_id,
            "residue_index": residue_index - 1,
            "atom_name": atom_name,
        }

    def _api_ligand_atom(self, chain_id: str, atom_name: str) -> ApiDict:
        """Convert a ligand atom reference, rejecting SMILES ligands."""
        seq = self._get_sequence_by_id(chain_id)
        if not isinstance(seq, Ligand):
            raise ValueError(
                f"Atom-level reference with atom name requires a ligand chain, got {chain_id!r}."
            )
        if seq.smiles:
            raise ValueError(
                f"Atom-level ligand references are not supported for SMILES ligand chain {chain_id!r}; use CCD ligands."
            )
        return {
            "type": "ligand_atom",
            "chain_id": chain_id,
            "atom_name": atom_name,
        }

    def _api_atom(self, atom: tuple[str, int, str]) -> ApiDict:
        """Convert a package bond atom to a Boltz API atom reference."""
        chain_id, residue_index, atom_name = atom
        seq = self._get_sequence_by_id(chain_id)
        if isinstance(seq, Ligand):
            return self._api_ligand_atom(chain_id, atom_name)
        return self._api_polymer_atom(chain_id, residue_index, atom_name)

    def _api_contact_token(self, token: tuple[str, int | str]) -> ApiDict:
        """Convert a package contact token to a Boltz API contact token."""
        chain_id, residue_or_atom = token
        seq = self._get_sequence_by_id(chain_id)
        if isinstance(residue_or_atom, int):
            if isinstance(seq, Ligand):
                raise ValueError(
                    f"Ligand contact token for chain {chain_id!r} must use an atom name, not a residue index."
                )
            if residue_or_atom < 1:
                raise ValueError(
                    f"Residue indices are 1-based and must be greater than zero: {residue_or_atom}"
                )
            return {
                "type": "polymer_contact",
                "chain_id": chain_id,
                "residue_index": residue_or_atom - 1,
            }

        if not isinstance(seq, Ligand):
            raise ValueError(
                f"Polymer contact token for chain {chain_id!r} must use a residue index, not an atom name."
            )
        if seq.smiles:
            raise ValueError(
                f"Atom-level ligand references are not supported for SMILES ligand chain {chain_id!r}; use CCD ligands."
            )
        return {
            "type": "ligand_contact",
            "chain_id": chain_id,
            "atom_name": residue_or_atom,
        }

    def _api_bond(self, bond: Bond) -> ApiDict:
        """Convert a bond constraint to the Boltz API top-level bonds list."""
        return {
            "atom1": self._api_atom(bond.atom1),
            "atom2": self._api_atom(bond.atom2),
        }

    def _api_contact(self, contact: Contact) -> ApiDict:
        """Convert a contact constraint to the Boltz API constraint shape."""
        d: ApiDict = {
            "type": "contact",
            "token1": self._api_contact_token(contact.token1),
            "token2": self._api_contact_token(contact.token2),
            "max_distance_angstrom": contact.max_distance,
        }
        if contact.force:
            d["force"] = contact.force
        return d

    def _api_pocket(self, pocket: Pocket) -> ApiDict:
        """Convert a pocket constraint to the Boltz API constraint shape."""
        contact_residues: dict[str, list[int]] = {}
        for contact_chain_id, residue_index in pocket.contacts:
            seq = self._get_sequence_by_id(contact_chain_id)
            if isinstance(seq, Ligand) or not isinstance(residue_index, int):
                raise ValueError(
                    "Pocket contacts for Boltz API output must be polymer residue positions."
                )
            if residue_index < 1:
                raise ValueError(
                    f"Residue indices are 1-based and must be greater than zero: {residue_index}"
                )
            contact_residues.setdefault(contact_chain_id, []).append(residue_index - 1)

        d: ApiDict = {
            "type": "pocket",
            "binder_chain_id": pocket.binder,
            "contact_residues": contact_residues,
            "max_distance_angstrom": pocket.max_distance,
        }
        if pocket.force:
            d["force"] = pocket.force
        return d

    def _api_template(self, template: Template) -> ApiDict:
        """Convert a local template definition to the Boltz API template shape."""
        if not template.chain_id:
            raise ValueError("Boltz API templates require at least one chain_id.")

        template_ids = template.template_id or template.chain_id
        if len(template.chain_id) != len(template_ids):
            raise ValueError(
                "Boltz API templates require matching chain_id and template_id lengths."
            )

        path_or_url = template.cif or template.pdb
        if path_or_url is None:
            raise ValueError("Template must define either a CIF or PDB source.")

        d: ApiDict = {
            "template_structure": self._file_source(
                path_or_url,
                {".cif": "chemical/x-cif", ".pdb": "chemical/x-pdb"},
            ),
            "template_chains": [
                {"input_chain_id": input_id, "template_chain_id": template_id}
                for input_id, template_id in zip(template.chain_id, template_ids)
            ],
        }
        if template.force:
            d["force_threshold_angstroms"] = template.threshold
        return d

    def _api_binding(self, protein_binder_chain_ids: str | list[str] | None) -> ApiDict | None:
        """Return the optional Boltz API binding object."""
        affinities = [prop for prop in self.properties if isinstance(prop, Affinity)]
        if protein_binder_chain_ids is not None and affinities:
            raise ValueError(
                "Cannot request both ligand-protein affinity and protein-protein binding in one Boltz API input."
            )

        if protein_binder_chain_ids is not None:
            binder_ids = (
                [protein_binder_chain_ids]
                if isinstance(protein_binder_chain_ids, str)
                else protein_binder_chain_ids
            )
            if not binder_ids:
                raise ValueError("protein_binder_chain_ids cannot be empty.")
            for chain_id in binder_ids:
                seq = self._get_sequence_by_id(chain_id)
                if not isinstance(seq, ProteinChain):
                    raise ValueError(
                        f"Protein-protein binding requires protein binder chains, got {chain_id!r}."
                    )
            return {
                "type": "protein_protein_binding",
                "binder_chain_ids": binder_ids,
            }

        if not affinities:
            return None

        binder = affinities[0].binder
        binder_seq = self._get_sequence_by_id(binder)
        if not isinstance(binder_seq, Ligand):
            raise ValueError(
                f"Ligand-protein binding requires a ligand binder chain, got {binder!r}."
            )
        if len(binder_seq.ids) != 1:
            raise ValueError(
                "Ligand-protein binding requires the ligand entity to have exactly one chain ID."
            )
        for seq in self.sequences:
            if not isinstance(seq, (ProteinChain, Ligand)):
                raise ValueError(
                    "Ligand-protein binding in the Boltz API only supports complexes containing proteins and ligands."
                )
        return {
            "type": "ligand_protein_binding",
            "binder_chain_id": binder,
        }

    @staticmethod
    def _validate_boltz_api_options(
        num_samples: int | None,
        recycling_steps: int | None,
        sampling_steps: int | None,
        step_scale: float | None,
    ) -> None:
        """Validate Boltz API prediction options."""
        if num_samples is not None and not (1 <= num_samples <= 10):
            raise ValueError("num_samples must be between 1 and 10.")
        if recycling_steps is not None and recycling_steps < 1:
            raise ValueError("recycling_steps must be greater than or equal to 1.")
        if sampling_steps is not None and sampling_steps < 50:
            raise ValueError("sampling_steps must be greater than or equal to 50.")
        if step_scale is not None and not (1.3 <= step_scale <= 2.0):
            raise ValueError("step_scale must be between 1.3 and 2.0.")

    def to_boltz_api_input(
        self,
        *,
        num_samples: int | None = None,
        recycling_steps: int | None = None,
        sampling_steps: int | None = None,
        step_scale: float | None = None,
        protein_binder_chain_ids: str | list[str] | None = None,
    ) -> dict[str, Any]:
        """Convert the job to a Boltz API structure-and-binding input object.

        The returned dictionary is the `input` body for
        `client.predictions.structure_and_binding.start(model=..., input=...)`.
        Unlike the YAML schema, the Boltz API uses 0-based residue indices.
        """
        self._check_ids()
        if not self.sequences:
            raise ValueError("Empty list of sequences.")

        self._validate_boltz_api_options(
            num_samples, recycling_steps, sampling_steps, step_scale
        )

        d: ApiDict = {
            "entities": [self._api_entity(seq) for seq in self.sequences],
        }

        if binding := self._api_binding(protein_binder_chain_ids):
            d["binding"] = binding

        bonds = [
            self._api_bond(constraint)
            for constraint in self.constraints
            if isinstance(constraint, Bond)
        ]
        if bonds:
            d["bonds"] = bonds

        constraints = [
            self._api_contact(constraint)
            if isinstance(constraint, Contact)
            else self._api_pocket(constraint)
            for constraint in self.constraints
            if isinstance(constraint, (Contact, Pocket))
        ]
        if constraints:
            d["constraints"] = constraints

        if self.templates:
            if len(self.templates) > 4:
                raise ValueError("Boltz API supports at most 4 templates.")
            d["templates"] = [self._api_template(template) for template in self.templates]

        model_options: ApiDict = {}
        if recycling_steps is not None:
            model_options["recycling_steps"] = recycling_steps
        if sampling_steps is not None:
            model_options["sampling_steps"] = sampling_steps
        if step_scale is not None:
            model_options["step_scale"] = step_scale
        if model_options:
            d["model_options"] = model_options

        if num_samples is not None:
            d["num_samples"] = num_samples

        return d

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

    def to_yaml(self) -> str:
        """Convert the job to a YAML string."""
        return yaml.dump(
            self.to_dict(), Dumper=IndentedDumper, indent=2, sort_keys=False
        )

    def write_yaml(self, filename: str) -> None:
        """Write the job to a YAML file as input for Boltz-2."""
        with open(filename, "w") as f:
            f.write(self.to_yaml())

    @classmethod
    def from_yaml(cls, filename: str, jobname: str = "boltz_job") -> Self:
        """Read a job from a YAML file."""
        sequence_match: dict[str, type[Sequence]] = {
            "protein": ProteinChain,
            "dna": DnaChain,
            "rna": RnaChain,
            "ligand": Ligand,
        }

        constraint_match: dict[str, type[Constraint]] = {
            "bond": Bond,
            "contact": Contact,
            "pocket": Pocket,
        }

        property_match: dict[str, type[Property]] = {
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

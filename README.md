# boltzjobs

`boltzjobs` is a Python package designed to streamline the process of creating YAML input files for [Boltz-2](https://github.com/jwohlwend/boltz).

It provides an intuitive, object-oriented interface for defining molecular sequences, constraints, and templates, and generates a correctly formatted YAML file.

## Features

- Define **protein**, **DNA**, and **RNA chains** with support for residue modifications.
- Add **ligands** using either CCD codes or SMILES strings.
- Specify geometric **constraints**, including bonds, contacts, and binding pockets.
- Define mmCIF files as **structural templates**.
- Request property calculations, such as **binding affinity**.

## Installation

```bash
pip install git+https://github.com/ugSUBMARINE/boltzjobs.git
```

## Usage

The `Job` class is the main entry point for building a Boltz-2 input file. You can create a job, add different molecular components, and then write the final configuration to a YAML file.

```python
from boltzjobs import Job

# 1. Create a new job
job = Job(name="Herceptin-HER2_complex_design")

# 2. Add molecular components
# Add a protein chain (HER2)
her2_chain = job.add_protein_chain(sequence="MELAALCRWGLLLALLPPGAASTQVCTGTDMKLRLPASPETHLDMLRHLYQGCQVVQGNLELTYLPTNASLSFLQDIQEVQGYVLIAHNQVRQVPLQRLRIVRGTQLFEDNYALAVLDNGDPLNNTTPVTGASPGGLRELQLRSLTEILKGGVLIQRNPQLCYQDTILWKDIFHKNNQLALTLIDTNRSR", ids="A")

# Add a second protein chain (Herceptin heavy chain)
herceptin_h = job.add_protein_chain(sequence="EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS", ids="H")

# Add a ligand using a SMILES string
# The library automatically handles single-quoting the SMILES string in the YAML output
job.add_ligand(smiles="CC(=O)N[C@@H](C)C(=O)O", ids="L")

# 3. Define constraints and properties
# Define a binding pocket for the ligand, specifying contact residues on chain A
pocket = job.add_pocket(binder="L", max_distance=6.0)
pocket.add_contact_token("A", 8)
pocket.add_contact_token("A", 12)

# Request an affinity calculation for the ligand
job.request_affinity(binder="L")

# 4. Add structural templates
job.add_template(cif="templates/1n8z.cif", chain_id="A")

# 5. Write the final YAML file
# The output will be correctly formatted for Boltz-2.
job.write_yaml("boltz_input.yaml")

print("Successfully generated boltz_input.yaml")
print("-" * 20)
print(job) # You can also print the job object for a summary
```

This will generate a `boltz_input.yaml` file with the following content:

```yaml
version: 1
sequences:
  - protein:
      id: A
      sequence: MELAALCRWGLLLALLPPGAASTQVCTGTDMKLRLPASPETHLDMLRHLYQGCQVVQGNLELTYLPTNASLSFLQDIQEVQGYVLIAHNQVRQVPLQRLRIVRGTQLFEDNYALAVLDNGDPLNNTTPVTGASPGGLRELQLRSLTEILKGGVLIQRNPQLCYQDTILWKDIFHKNNQLALTLIDTNRSR
  - protein:
      id: H
      sequence: EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS
  - ligand:
      id: L
      smiles: 'CC(=O)N[C@@H](C)C(=O)O'
constraints:
  - pocket:
      binder: L
      contacts: [[A, 8], [A, 12]]
      max_distance: 6.0
templates:
  - cif: templates/1n8z.cif
    chain_id: A
properties:
  - affinity:
      binder: L
```

## License

This project is licensed under the MIT License.

## References

- [Boltz-2 GitHub Repository](https://github.com/jwohlwend/boltz)
- [Boltz-2 Input File Schema Documentation](https://github.com/jwohlwend/boltz/blob/main/docs/prediction.md)
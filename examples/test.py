# -*- coding: utf-8 -*-
from boltzjobs import Job

# 1. Create a new job
job = Job(name="Herceptin-HER2_complex_design")

# 2. Add molecular components
# Add a protein chain (HER2)
her2_chain = job.add_protein_chain(
    sequence="MELAALCRWGLLLALLPPGAASTQVCTGTDMKLRLPASPETHLDMLRHLYQGCQVVQGNLELTYLPTNASLSFLQDIQEVQGYVLIAHNQVRQVPLQRLRIVRGTQLFEDNYALAVLDNGDPLNNTTPVTGASPGGLRELQLRSLTEILKGGVLIQRNPQLCYQDTILWKDIFHKNNQLALTLIDTNRSR",
    ids="A",
)

# Add a second protein chain (Herceptin heavy chain)
herceptin_h = job.add_protein_chain(
    sequence="EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS",
    ids="H",
)

# Add a ligand using a SMILES string
# The library automatically handles single-quoting the SMILES string in the YAML output
_ = job.add_ligand(smiles="CC(=O)N[C@@H](C)C(=O)O", ids="L")

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
print(job)  # You can also print the job object for a summary


# from boltzjobs import Job

# job = Job.from_yaml("examples/prot_no_msa.yaml")
# print(job)
# job.write_yaml(f"{job.name}.yaml")

# job = Job("boltz_job")
# p = job.add_protein_chain("ACDEFGHIKLMNPQRSTVWY", count=2, cyclic=True, msa="test.a3m")
# p.add_modification("XXX", 2)
# p.add_modification("YYY", 3)
#
# job.add_ligand(ccd="HEM", ids="X")
# job.add_ligand(smiles="N[C@@H](Cc1ccc(O)cc1)C(=O)O")
# job.request_affinity("C")
#
# job.add_bond("A", 5, "OE2", "X", 1, "CB")
# pocket = job.add_pocket("C", contact_tokens=[("A", 5)], max_distance=5.0)
# pocket.add_contact_token("A", 1)
# pocket.add_contact_token("B", 2)
#
# job.add_contact("C", 1, "A", 1, max_distance=4.0)
#
# job.add_template("test.cif", "A")
# job.add_template("test_2.cif", "B", "C")
# job.add_template("test_3.cif")
#
# d = job.to_dict()
# print(d)
#
# job.write_yaml(f"{job.name}.yaml")
#
# print(job)

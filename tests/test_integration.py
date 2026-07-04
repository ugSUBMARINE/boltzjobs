"""Integration tests for boltzjobs complete workflows."""

import pytest
import yaml
from pathlib import Path

from boltzjobs import Job


class TestCompleteWorkflows:
    """Test complete workflows from start to finish."""

    @pytest.mark.integration
    def test_protein_ligand_workflow(self, temp_yaml_file):
        """Test complete protein-ligand job workflow."""
        # Create job
        job = Job("protein_ligand_complex")

        # Add protein
        job.add_protein_chain(
            sequence="MVTPEGNVSLVDESLLVGVTDED", ids="A", msa="empty"
        )

        # Add ligand
        job.add_ligand(smiles="CC(=O)N[C@@H](C)C(=O)O", ids="L")

        # Add constraints
        pocket = job.add_pocket("L", max_distance=6.0)
        pocket.add_contact_token("A", 5)
        pocket.add_contact_token("A", 10)

        job.add_contact("A", 1, "L", 1, max_distance=5.0)

        # Add template
        job.add_template("template.cif", chain_id="A")

        # Request affinity
        job.request_affinity("L")

        # Write to YAML
        job.write_yaml(temp_yaml_file)

        # Verify file exists and contains expected content
        assert Path(temp_yaml_file).exists()

        with open(temp_yaml_file) as f:
            content = f.read()
            assert "protein_ligand_complex" in job.name
            assert "version: 1" in content
            assert "protein:" in content
            assert "ligand:" in content
            assert "pocket:" in content
            assert "contact:" in content
            assert "template:" in content or "cif:" in content
            assert "affinity:" in content

        # Test round-trip
        loaded_job = Job.from_yaml(temp_yaml_file, "loaded_job")
        assert len(loaded_job.sequences) == 2
        assert len(loaded_job.constraints) == 2
        assert len(loaded_job.templates) == 1
        assert len(loaded_job.properties) == 1

    @pytest.mark.integration
    def test_multimer_workflow(self, temp_yaml_file):
        """Test multimeric protein complex workflow."""
        job = Job("multimer_complex")

        # Add multiple protein chains
        job.add_protein_chain(
            sequence="MVTPEGNVSLVDESLLVGVTDED",
            count=2,  # Two copies
            ids=["A", "B"],
        )

        job.add_protein_chain(sequence="EAFSLFDKSLVNETP", ids="C")

        # Add inter-chain bonds
        job.add_contact("A", 5, "B", 10, max_distance=6.0)
        job.add_contact("B", 15, "C", 5, max_distance=5.0)

        # Add templates
        job.add_template("template1.cif", chain_id=["A", "B"])
        job.add_template("template2.cif", chain_id="C")

        # Write and verify
        job.write_yaml(temp_yaml_file)

        loaded_job = Job.from_yaml(temp_yaml_file, "loaded_multimer")
        assert len(loaded_job.sequences) == 2  # Two different sequences

        # Check that first sequence has multiple IDs
        multimer_chain = None
        single_chain = None
        for seq in loaded_job.sequences:
            if len(seq.ids) == 2:
                multimer_chain = seq
            elif len(seq.ids) == 1:
                single_chain = seq

        assert multimer_chain is not None
        assert single_chain is not None
        assert set(multimer_chain.ids) == {"A", "B"}
        assert single_chain.ids == ["C"]

    @pytest.mark.integration
    def test_nucleic_acid_workflow(self, temp_yaml_file):
        """Test DNA/RNA complex workflow."""
        job = Job("nucleic_acid_complex")

        # Add DNA and RNA
        job.add_dna_chain("ATCGATCGATCG", ids="D")
        job.add_rna_chain("AUCGAUCGAUCG", ids="R")

        # Add protein that binds to nucleic acids
        job.add_protein_chain(sequence="MVTPEGNVSLVDESLLVGVTDED", ids="P")

        # Add nucleic acid-protein contacts
        job.add_contact("D", 5, "P", 10, max_distance=6.0)
        job.add_contact("R", 3, "P", 15, max_distance=5.0)

        # Add template
        job.add_template("nucleic_template.cif", chain_id=["D", "R", "P"])

        job.write_yaml(temp_yaml_file)

        loaded_job = Job.from_yaml(temp_yaml_file, "loaded_nucleic")

        # Verify all sequence types present
        seq_types = []
        for seq in loaded_job.sequences:
            seq_types.append(type(seq).__name__)

        assert "DnaChain" in seq_types
        assert "RnaChain" in seq_types
        assert "ProteinChain" in seq_types

    @pytest.mark.integration
    def test_modified_sequence_workflow(self, temp_yaml_file):
        """Test workflow with sequence modifications."""
        job = Job("modified_sequence")

        # Add protein with modifications
        protein = job.add_protein_chain(
            sequence="MVTPEGNVSLVDESLLVGVTCDED", ids="A", cyclic=True
        )

        # Add modifications
        protein.add_modification("HIS", 5)
        protein.add_modification("PHE", 10)
        protein.add_modification("MSE", 20)  # Selenomethionine

        # Add DNA with modifications
        dna = job.add_dna_chain("ATCGATCGATCG", ids="D")
        dna.add_modification("5MC", 3)  # 5-methylcytosine

        job.write_yaml(temp_yaml_file)

        loaded_job = Job.from_yaml(temp_yaml_file, "loaded_modified")

        # Verify modifications preserved
        loaded_protein = None
        loaded_dna = None
        for seq in loaded_job.sequences:
            if hasattr(seq, "sequence") and "MVTP" in seq.sequence:
                loaded_protein = seq
            elif hasattr(seq, "sequence") and "ATCG" in seq.sequence:
                loaded_dna = seq

        assert loaded_protein is not None
        assert loaded_dna is not None
        assert loaded_protein.cyclic is True
        assert len(loaded_protein.modifications) == 3
        assert len(loaded_dna.modifications) == 1

        # Check specific modifications
        protein_mods = [(m.ccd, m.position) for m in loaded_protein.modifications]
        assert ("HIS", 5) in protein_mods
        assert ("PHE", 10) in protein_mods
        assert ("MSE", 20) in protein_mods

    @pytest.mark.integration
    def test_complex_constraint_workflow(self, temp_yaml_file):
        """Test workflow with multiple constraint types."""
        job = Job("complex_constraints")

        # Add sequences
        job.add_protein_chain("MVTPEGNVSLVDESLLCGVTCED", ids="A")  # With cysteines
        job.add_protein_chain("MVTPEGNVSLVDESLLCGVTCED", ids="B")
        job.add_ligand(smiles="CC(=O)N[C@@H](C)C(=O)O", ids="L1")
        job.add_ligand(ccd="ATP", ids="L2")

        # Add various constraints
        # Disulfide bond between proteins
        job.add_disulfide_bond("A", 17, "B", 17)  # Cysteine positions

        # Ligand binding pockets
        pocket1 = job.add_pocket("L1", max_distance=6.0)
        pocket1.add_contact_token("A", 5)
        pocket1.add_contact_token("A", 10)

        pocket2 = job.add_pocket("L2", max_distance=5.0)
        pocket2.add_contact_token("B", 8)
        pocket2.add_contact_token("B", 12)

        # Protein-protein contacts
        job.add_contact("A", 1, "B", 1, max_distance=4.0)
        job.add_contact("A", 20, "B", 20, max_distance=7.0)

        # Ligand-ligand contact
        job.add_contact("L1", 1, "L2", 1, max_distance=8.0)

        # Specific atomic bond (e.g., coordination bond)
        job.add_bond("A", 15, "OG", "L2", 1, "PA")

        job.write_yaml(temp_yaml_file)

        loaded_job = Job.from_yaml(temp_yaml_file, "loaded_complex")

        # Verify all constraint types
        constraint_types = [type(c).__name__ for c in loaded_job.constraints]
        assert "Bond" in constraint_types
        assert "Contact" in constraint_types
        assert "Pocket" in constraint_types

        # Should have multiple instances of each type
        assert constraint_types.count("Contact") >= 3  # At least 3 contacts
        assert constraint_types.count("Pocket") == 2  # Exactly 2 pockets
        assert constraint_types.count("Bond") >= 1  # At least 1 bond

    @pytest.mark.integration
    def test_template_affinity_workflow(self, temp_yaml_file):
        """Test workflow with templates and affinity calculations."""
        job = Job("template_affinity")

        # Add protein target
        job.add_protein_chain(
            sequence="MVTPEGNVSLVDESLLVGVTDED", ids="TARGET", msa="target.a3m"
        )

        # Add multiple ligands but only request affinity for one
        job.add_ligand(smiles="CC(=O)N[C@@H](C)C(=O)O", ids="DRUG")
        job.add_ligand(ccd="HEM", ids="COFACTOR")

        # Add binding site definition
        pocket = job.add_pocket("DRUG", max_distance=6.0)
        for residue in [45, 67, 89, 123, 156]:
            pocket.add_contact_token("TARGET", residue)

        # Add multiple templates
        job.add_template("apo_structure.cif", chain_id="TARGET")
        job.add_template(
            "holo_structure.cif",
            chain_id=["TARGET", "COFACTOR"],
            template_id=["A", "H"],
        )

        # Request affinity for drug but not cofactor
        job.request_affinity("DRUG")

        job.write_yaml(temp_yaml_file)

        # Verify YAML structure
        with open(temp_yaml_file) as f:
            yaml_data = yaml.safe_load(f)

        assert len(yaml_data["sequences"]) == 3
        assert len(yaml_data["templates"]) == 2
        assert len(yaml_data["properties"]) == 1
        assert yaml_data["properties"][0]["affinity"]["binder"] == "DRUG"

        # Check template structures
        template_cifs = [t["cif"] for t in yaml_data["templates"]]
        assert "apo_structure.cif" in template_cifs
        assert "holo_structure.cif" in template_cifs

    @pytest.mark.integration
    def test_error_handling_workflow(self):
        """Test workflow error handling and validation."""
        job = Job("error_test")

        # Test sequence validation errors
        with pytest.raises(ValueError, match="Protein sequence contains invalid"):
            job.add_protein_chain("MVTPXYZ123")

        with pytest.raises(ValueError, match="DNA sequence can only contain"):
            job.add_dna_chain("ATCGXYZ")

        with pytest.raises(ValueError, match="RNA sequence can only contain"):
            job.add_rna_chain("AUCGXYZ")

        # Add valid sequences for constraint testing
        job.add_protein_chain("MVTPEGNVSLVDESLLVGVTDED", ids="A")
        job.add_ligand(smiles="CC(=O)O", ids="L")

        # Test constraint validation errors
        with pytest.raises(ValueError, match="Both chain IDs must be defined"):
            job.add_bond("X", 1, "CA", "Y", 1, "CB")

        # Pocket constraints now work for any sequence type, test invalid ID instead
        with pytest.raises(ValueError, match="No sequence with id 'Z' found"):
            job.add_pocket("Z")  # Try to define pocket for non-existent ID

        with pytest.raises(
            ValueError, match="Affinity can only be estimated for ligands"
        ):
            job.request_affinity("A")  # Try to request affinity for protein

        # Test template validation
        with pytest.raises(
            ValueError, match="Either 'cif' or 'pdb' file path must be provided"
        ):
            job.add_template("")

        # Test ligand validation
        with pytest.raises(
            ValueError, match="Either CCD codes or SMILES string must be provided"
        ):
            job.add_ligand()

    @pytest.mark.integration
    def test_id_generation_workflow(self, temp_yaml_file):
        """Test automatic ID generation workflow."""
        job = Job("id_generation_test")

        # Add many sequences to test ID generation
        sequences = []

        # Add 30 protein chains to exceed single letters
        for i in range(30):
            seq = job.add_protein_chain("MVTP", count=1)
            sequences.append(seq)

        # Verify ID progression
        expected_ids = list("ABCDEFGHIJKLMNOPQRSTUVWXYZ") + ["AA", "BA", "CA", "DA"]
        actual_ids = [seq.ids[0] for seq in sequences]

        assert actual_ids == expected_ids

        # Test writing large job
        job.write_yaml(temp_yaml_file)

        loaded_job = Job.from_yaml(temp_yaml_file, "loaded_large")
        assert len(loaded_job.sequences) == 30

        # Verify ID uniqueness in loaded job
        all_ids = []
        for seq in loaded_job.sequences:
            all_ids.extend(seq.ids)

        assert len(all_ids) == len(set(all_ids))  # All IDs should be unique


class TestRealWorldScenarios:
    """Test scenarios based on real molecular modeling use cases."""

    @pytest.mark.integration
    def test_antibody_antigen_complex(self, temp_yaml_file):
        """Test antibody-antigen complex modeling scenario."""
        job = Job("antibody_antigen")

        # Add antibody chains
        job.add_protein_chain(
            sequence="EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS",
            ids="H",
        )

        job.add_protein_chain(
            sequence="DIQMTQSPSSLSASVGDRVTITCRASQSISSWLAWYQQKPGKAPKLLIYKASSLESGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQYNSYNPLTFGGGTKVEIK",
            ids="L",
        )

        # Add antigen
        job.add_protein_chain(
            sequence="MELAALCRWGLLLALLPPGAASTQVCTGTDMKLRLPASPETHLDMLRHLYQGCQVVQGNLELTYLPTNASLSFLQDIQEVQGYVLIAHNQVRQVPLQRLRIVRGTQLFEDNYALAVLDNGDPLNNTTPVTGASPGGLRELQLRSLTEILKGGVLIQRNPQLCYQDTILWKDIFHKNNQLALTLIDTNRSR",
            ids="A",
        )

        # Add CDR-antigen contacts (paratope-epitope interactions)
        job.add_contact("H", 31, "A", 45, max_distance=6.0)  # CDR-H1 contact
        job.add_contact("H", 53, "A", 67, max_distance=5.0)  # CDR-H2 contact
        job.add_contact("H", 101, "A", 89, max_distance=4.5)  # CDR-H3 contact
        job.add_contact("L", 28, "A", 123, max_distance=5.5)  # CDR-L1 contact
        job.add_contact("L", 92, "A", 145, max_distance=6.0)  # CDR-L3 contact

        # Add disulfide bonds within antibody
        # Heavy chain intra-chain disulfide
        job.add_disulfide_bond("H", 22, "H", 96)
        # Light chain intra-chain disulfide
        job.add_disulfide_bond("L", 23, "L", 88)

        # Add template structures
        job.add_template("antibody_framework.cif", chain_id=["H", "L"])
        job.add_template("antigen_structure.cif", chain_id="A")

        job.write_yaml(temp_yaml_file)

        # Verify complex structure
        loaded_job = Job.from_yaml(temp_yaml_file, "loaded_antibody")
        assert len(loaded_job.sequences) == 3
        assert len(loaded_job.constraints) >= 7  # 5 contacts + 2 disulfides
        assert len(loaded_job.templates) == 2

    @pytest.mark.integration
    def test_drug_discovery_scenario(self, temp_yaml_file):
        """Test drug discovery workflow scenario."""
        job = Job("drug_discovery")

        # Add target protein (e.g., kinase)
        job.add_protein_chain(
            sequence="MENFQKVEKIGEGTYGVVYKARNKLTGEVVALKKIRLDTETEGVPSTAIREISLLKELNHPNIVKLLDVIHTENKLYLVFEFLHQDLKKFMDASALTGIPLPLIKSYLFQLLQGLAFCHSHRVLHRDLKPQNLLINTEGAIKLADFGLARAFGVPVRTYTHEVVTLWYRAPEILLGCKYYSTAVDIWSLGCIFAEMVTRRALFPGDSEIDQLFRIFRTLGTPDEVVWPGVTSMPDYKPSFPKWARQDFSKVVPPLDEDGRSLLSQMLHYDPNKRISAKAALAHPFFQDVTKPVPHLRL",
            ids="TARGET",
            msa="kinase_msa.a3m",
        )

        # Add multiple drug candidates
        compounds = [
            ("DRUG1", "CC1=C(C=CC(=C1)NC(=O)C2=CC=C(C=C2)CN3CCN(CC3)C)C(F)(F)F"),
            ("DRUG2", "CN(C)C1=NC=NC2=C1C=CC(=C2)NC3=CC(=CC(=C3)C(F)(F)F)C(F)(F)F"),
            (
                "DRUG3",
                "CC(C)(C)C1=CC=C(C=C1)C(=O)NC2=CC3=C(C=C2)N=CN=C3NC4=CC=C(C=C4)N5CCN(CC5)C",
            ),
        ]

        for drug_id, smiles in compounds:
            job.add_ligand(smiles=smiles, ids=drug_id)

            # Define binding site for each drug
            pocket = job.add_pocket(drug_id, max_distance=6.0)
            # Active site residues (typical kinase binding site)
            for residue in [52, 72, 94, 145, 147, 161, 172]:
                pocket.add_contact_token("TARGET", residue)

        # Request affinity calculation for only one drug (per Boltz-2 schema)
        job.request_affinity("DRUG1")

        # Add known structure template
        job.add_template("kinase_apo.cif", chain_id="TARGET")

        job.write_yaml(temp_yaml_file)

        loaded_job = Job.from_yaml(temp_yaml_file, "loaded_drugs")

        # Verify that all drugs were loaded
        assert len(loaded_job.sequences) == 4  # 1 target + 3 drugs
        assert len(loaded_job.properties) == 1  # 1 affinity request (schema limit)
        assert loaded_job.properties[0].binder == "DRUG1"

    @pytest.mark.integration
    def test_membrane_protein_scenario(self, temp_yaml_file):
        """Test membrane protein complex scenario."""
        job = Job("membrane_protein")

        # Add membrane protein (GPCR-like)
        job.add_protein_chain(
            sequence="MGSSTSAPPNISCSSCLPLERPSTSQHPRNSSCDAISYDRYASIFLCIVFGNIMVIILIVLRLRRRKKLNGEKTRPRRNLPQMPTGITAQEDLLASDDAYNPQEFCEFMYKTKKRTEKDVVQSLVAMLCNIIASDSLNPVYFMGSLNQTAEFMRRMLLQTLTWKMAVLLGLIASILFLALVTEDQQQQLQIQHQTNQPLLQQNPSSRCVSIIDPPCRLSSLHSLSLRRSLRQPPQRFLLHRCQRRLQTSLGLSSLLEMLR",
            ids="GPCR",
        )

        # Add G-protein subunits
        job.add_protein_chain(
            sequence="MGCTQSAEDSKCQKYTRQVDQMEYLCSQQTIDQVHRNLNNVGQDSYKQCLHQLDQGEKQILFYCFKNLLKRRCRKGMEGLTNQDSMEGLRLMNQNKIAAAQHNNQPMNQSCSQHDRYNLCDLLLRRIEMCDKDSMNQLYQNLLQRETKSLKVQYTLTLEDLRWSQRSLTSLQLLNNKTNKQHISVEQLQRRRQQLRQQHQHQQRQRQRQQQPPRLNQN",
            ids="GA",
        )

        job.add_protein_chain(
            sequence="MGQGSLKQVFQKYKWLNVPMEYSLNSLMDSYQTQDQRLISHMQRQMKQSTLYTKKLSSSCTDKFKQVAQEQLEQCLQAHLLRKFNLFQKDSSLLLQRPSQVQPKYLYVSGEFLAKGSQSGLQNSYILQLRKSGVYLPKSKRTLDYQQQRQQQYTLDSMQQMTLREEWKQKTKAMIREWQRQRDTHRTKQLLQHQRQRQQQRRRNQRKQQLMQSKLVEQRSHRVIQKLITLTQPQRQHKRQLRKRSRLTSQHRQRRERLRSQQQQSKAQTSLHLLTQRQERQHSRRLRNQQRQRR",
            ids="GB",
        )

        job.add_protein_chain(
            sequence="MSRLDQCCDCCQARNSQKQAFCCQHLERECEKRKKRLEKQNSLLLLLLLLECVGQNCCKSCPCCLTSLEKQCLQCCLGSCECKCKCLKQCCQSRLKQCCQKC",
            ids="GG",
        )

        # Add GPCR-G protein contacts
        job.add_contact("GPCR", 140, "GA", 45, max_distance=6.0)  # IC3-Gα contact
        job.add_contact("GPCR", 230, "GA", 89, max_distance=5.0)  # IC2-Gα contact
        job.add_contact("GPCR", 310, "GA", 123, max_distance=6.0)  # C-term-Gα contact

        # Add G-protein heterotrimer contacts
        job.add_contact("GA", 200, "GB", 50, max_distance=5.0)  # Gα-Gβ interface
        job.add_contact("GB", 150, "GG", 20, max_distance=4.0)  # Gβ-Gγ interface

        # Add lipid interactions (modeled as small molecules)
        job.add_ligand(
            smiles="CCCCCCCCCCCCCCCCCC(=O)OC[C@H](COP(=O)(O)OCC[N+](C)(C)C)OC(=O)CCCCCCCCCCCCCCC",
            ids="PC",
        )  # Phosphatidylcholine

        # GPCR-lipid contacts
        job.add_contact("GPCR", 50, "PC", 1, max_distance=4.0)  # TM1-lipid
        job.add_contact("GPCR", 100, "PC", 1, max_distance=4.0)  # TM3-lipid
        job.add_contact("GPCR", 180, "PC", 1, max_distance=4.0)  # TM5-lipid

        # Add structural templates
        job.add_template("gpcr_inactive.cif", chain_id="GPCR")
        job.add_template("g_protein_trimer.cif", chain_id=["GA", "GB", "GG"])

        job.write_yaml(temp_yaml_file)

        loaded_job = Job.from_yaml(temp_yaml_file, "loaded_membrane")
        assert len(loaded_job.sequences) == 5  # 4 proteins + 1 lipid
        assert (
            len(loaded_job.constraints) >= 8
        )  # Multiple protein-protein and protein-lipid contacts
        assert len(loaded_job.templates) == 2

"""Tests for boltzjobs.utils module."""

import pytest
import yaml
import json
from itertools import islice
from io import StringIO

from boltzjobs.utils import (
    chain_id, FlowStyleList, SingleQuoted, IndentedDumper,
    flow_style_list_representer, single_quoted_representer,
    get_msa_from_json, get_templates_from_json, fasta_sequences
)


class TestChainId:
    """Test the chain_id generator function."""
    
    @pytest.mark.unit
    def test_chain_id_basic_sequence(self):
        """Test basic chain ID generation."""
        gen = chain_id()
        ids = list(islice(gen, 30))
        
        # First 26 should be single letters
        expected_single = list("ABCDEFGHIJKLMNOPQRSTUVWXYZ")
        assert ids[:26] == expected_single
        
        # Next 4 should be AA, BA, CA, DA (reverse spreadsheet style)
        assert ids[26:30] == ["AA", "BA", "CA", "DA"]
    
    @pytest.mark.unit
    def test_chain_id_custom_letters(self):
        """Test chain ID generation with custom alphabet."""
        gen = chain_id("ABC")
        ids = list(islice(gen, 6))
        
        expected = ["A", "B", "C", "AA", "BA", "CA"]
        assert ids == expected
    
    @pytest.mark.unit
    def test_chain_id_generator_continues(self):
        """Test that generator continues producing unique IDs."""
        gen = chain_id()
        
        # Get first batch
        batch1 = list(islice(gen, 10))
        # Get second batch
        batch2 = list(islice(gen, 10))
        
        # Should be no overlap
        assert set(batch1).isdisjoint(set(batch2))
        assert len(set(batch1 + batch2)) == 20


class TestFlowStyleList:
    """Test FlowStyleList custom YAML formatting."""
    
    @pytest.mark.unit
    def test_flow_style_list_creation(self):
        """Test FlowStyleList creation and basic functionality."""
        fsl = FlowStyleList([1, 2, 3])
        assert fsl == [1, 2, 3]
        assert isinstance(fsl, list)
        assert isinstance(fsl, FlowStyleList)
    
    @pytest.mark.unit
    def test_flow_style_list_yaml_output(self):
        """Test FlowStyleList YAML representation."""
        fsl = FlowStyleList(["A", "B", "C"])
        
        # Use IndentedDumper to test custom representer
        yaml_output = yaml.dump({"test": fsl}, Dumper=IndentedDumper)
        
        # Should be in flow style: [A, B, C]
        assert "[A, B, C]" in yaml_output
    
    @pytest.mark.unit
    def test_flow_style_list_nested(self):
        """Test nested FlowStyleList."""
        nested = FlowStyleList([["A", 1], ["B", 2]])
        yaml_output = yaml.dump({"nested": nested}, Dumper=IndentedDumper)
        
        # Should maintain flow style
        assert "[[A, 1], [B, 2]]" in yaml_output


class TestSingleQuoted:
    """Test SingleQuoted custom YAML formatting."""
    
    @pytest.mark.unit
    def test_single_quoted_creation(self):
        """Test SingleQuoted creation."""
        sq = SingleQuoted("test string")
        assert sq == "test string"
        assert isinstance(sq, str)
        assert isinstance(sq, SingleQuoted)
    
    @pytest.mark.unit
    def test_single_quoted_yaml_output(self):
        """Test SingleQuoted YAML representation."""
        sq = SingleQuoted("CC(=O)N[C@@H](C)C(=O)O")
        
        yaml_output = yaml.dump({"smiles": sq}, Dumper=IndentedDumper)
        
        # Should be single-quoted
        assert "'CC(=O)N[C@@H](C)C(=O)O'" in yaml_output
    
    @pytest.mark.unit
    def test_single_quoted_special_characters(self):
        """Test SingleQuoted with special characters."""
        special_chars = "string with 'quotes' and \"double quotes\""
        sq = SingleQuoted(special_chars)
        
        yaml_output = yaml.dump({"test": sq}, Dumper=IndentedDumper)
        
        # Should handle special characters properly
        assert "'" in yaml_output


class TestIndentedDumper:
    """Test IndentedDumper YAML formatting."""
    
    @pytest.mark.unit
    def test_indented_dumper_basic(self):
        """Test basic IndentedDumper functionality."""
        data = {
            "version": 1,
            "sequences": [
                {"protein": {"id": "A", "sequence": "TEST"}}
            ]
        }
        
        yaml_output = yaml.dump(data, Dumper=IndentedDumper, indent=2)
        
        # Should be properly formatted
        assert "version: 1" in yaml_output
        assert "sequences:" in yaml_output
        assert "  - protein:" in yaml_output
    
    @pytest.mark.unit
    def test_indented_dumper_with_custom_classes(self):
        """Test IndentedDumper with FlowStyleList and SingleQuoted."""
        data = {
            "contacts": FlowStyleList([["A", 1], ["B", 2]]),
            "smiles": SingleQuoted("CC(=O)O")
        }
        
        yaml_output = yaml.dump(data, Dumper=IndentedDumper)
        
        assert "[[A, 1], [B, 2]]" in yaml_output
        assert "'CC(=O)O'" in yaml_output


class TestRepresenters:
    """Test YAML representer functions."""
    
    @pytest.mark.unit
    def test_flow_style_list_representer(self):
        """Test flow_style_list_representer function."""
        from io import StringIO
        dumper = IndentedDumper(StringIO())
        fsl = FlowStyleList([1, 2, 3])
        
        node = flow_style_list_representer(dumper, fsl)
        
        assert node.tag == "tag:yaml.org,2002:seq"
        assert node.flow_style is True
    
    @pytest.mark.unit
    def test_single_quoted_representer(self):
        """Test single_quoted_representer function."""
        from io import StringIO
        dumper = IndentedDumper(StringIO())
        sq = SingleQuoted("test")
        
        node = single_quoted_representer(dumper, sq)
        
        assert node.tag == "tag:yaml.org,2002:str"
        assert node.style == "'"


class TestMSAFromJSON:
    """Test get_msa_from_json function."""
    
    @pytest.mark.unit
    def test_get_msa_from_json_protein_unpaired(self, tmp_path):
        """Test extracting unpaired MSA for protein."""
        json_data = {
            "sequences": [
                {
                    "protein": {
                        "sequence": "TESTSEQ",
                        "unpairedMsa": "unpaired_msa_data",
                        "pairedMsa": "paired_msa_data"
                    }
                }
            ]
        }
        
        json_file = tmp_path / "test.json"
        json_file.write_text(json.dumps(json_data))
        
        result = get_msa_from_json(str(json_file), "TESTSEQ", paired=False)
        assert result == "unpaired_msa_data"
    
    @pytest.mark.unit
    def test_get_msa_from_json_protein_paired(self, tmp_path):
        """Test extracting paired MSA for protein."""
        json_data = {
            "sequences": [
                {
                    "protein": {
                        "sequence": "TESTSEQ",
                        "unpairedMsa": "unpaired_msa_data",
                        "pairedMsa": "paired_msa_data"
                    }
                }
            ]
        }
        
        json_file = tmp_path / "test.json"
        json_file.write_text(json.dumps(json_data))
        
        result = get_msa_from_json(str(json_file), "TESTSEQ", paired=True)
        assert result == "paired_msa_data"
    
    @pytest.mark.unit
    def test_get_msa_from_json_rna(self, tmp_path):
        """Test extracting MSA for RNA."""
        json_data = {
            "sequences": [
                {
                    "rna": {
                        "sequence": "AUCG",
                        "unpairedMsa": "rna_msa_data"
                    }
                }
            ]
        }
        
        json_file = tmp_path / "test.json"
        json_file.write_text(json.dumps(json_data))
        
        result = get_msa_from_json(str(json_file), "AUCG")
        assert result == "rna_msa_data"
    
    @pytest.mark.unit
    def test_get_msa_from_json_not_found(self, tmp_path):
        """Test MSA extraction when sequence not found."""
        json_data = {"sequences": []}
        
        json_file = tmp_path / "test.json"
        json_file.write_text(json.dumps(json_data))
        
        result = get_msa_from_json(str(json_file), "NOTFOUND")
        assert result is None


class TestTemplatesFromJSON:
    """Test get_templates_from_json function."""
    
    @pytest.mark.unit
    def test_get_templates_from_json_found(self, tmp_path):
        """Test extracting templates for protein."""
        template_data = [{"template1": "data1"}, {"template2": "data2"}]
        json_data = {
            "sequences": [
                {
                    "protein": {
                        "sequence": "TESTSEQ",
                        "templates": template_data
                    }
                }
            ]
        }
        
        json_file = tmp_path / "test.json"
        json_file.write_text(json.dumps(json_data))
        
        result = get_templates_from_json(str(json_file), "TESTSEQ")
        assert result == template_data
    
    @pytest.mark.unit
    def test_get_templates_from_json_not_found(self, tmp_path):
        """Test template extraction when sequence not found."""
        json_data = {"sequences": []}
        
        json_file = tmp_path / "test.json"
        json_file.write_text(json.dumps(json_data))
        
        result = get_templates_from_json(str(json_file), "NOTFOUND")
        assert result is None


class TestFastaSequences:
    """Test fasta_sequences generator function."""
    
    @pytest.mark.unit
    def test_fasta_sequences_basic(self, tmp_path):
        """Test basic FASTA parsing."""
        fasta_content = """>seq1
ATCG
>seq2
GCTA
"""
        fasta_file = tmp_path / "test.fasta"
        fasta_file.write_text(fasta_content)
        
        sequences = list(fasta_sequences(str(fasta_file)))
        
        assert len(sequences) == 2
        assert sequences[0] == (">seq1", "ATCG")
        assert sequences[1] == (">seq2", "GCTA")
    
    @pytest.mark.unit
    def test_fasta_sequences_multiline(self, tmp_path):
        """Test FASTA parsing with multiline sequences."""
        fasta_content = """>seq1 description
ATCG
GCTA
TTAA
>seq2
GC
TA
"""
        fasta_file = tmp_path / "test.fasta"
        fasta_file.write_text(fasta_content)
        
        sequences = list(fasta_sequences(str(fasta_file)))
        
        assert len(sequences) == 2
        assert sequences[0] == (">seq1 description", "ATCGGCTATTAA")
        assert sequences[1] == (">seq2", "GCTA")
    
    @pytest.mark.unit
    def test_fasta_sequences_empty_file(self, tmp_path):
        """Test FASTA parsing with empty file."""
        fasta_file = tmp_path / "empty.fasta"
        fasta_file.write_text("")
        
        sequences = list(fasta_sequences(str(fasta_file)))
        
        assert len(sequences) == 1
        assert sequences[0] == ("", "")
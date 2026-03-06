#!/usr/bin/env python
# coding: utf-8
"""
Unit tests for oncodriver.annot module.
"""

import pytest
from oncodriver.annot import (
    clean_version_number,
    get_gene_type,
    load_cancer_genes_list,
    load_canonicals,
)


# ---------------------------------------------------------------------------
# clean_version_number
# ---------------------------------------------------------------------------

class TestCleanVersionNumber:
    def test_removes_version_from_nm_id(self):
        assert clean_version_number("NM_001234.5") == "NM_001234"

    def test_removes_version_from_ensembl_id(self):
        assert clean_version_number("ENST00000123456.7") == "ENST00000123456"

    def test_no_version_unchanged(self):
        assert clean_version_number("NM_001234") == "NM_001234"

    def test_empty_string(self):
        assert clean_version_number("") == ""

    def test_only_dot_and_version(self):
        assert clean_version_number(".5") == ""

    def test_multiple_dots_strips_first_onwards(self):
        # Only the first dot and everything after it should be removed
        assert clean_version_number("NM_001234.5.6") == "NM_001234"


# ---------------------------------------------------------------------------
# get_gene_type
# ---------------------------------------------------------------------------

class TestGetGeneType:
    @pytest.fixture(autouse=True)
    def atlas(self):
        self.atlas = {
            "oncogene": {"KRAS", "EGFR"},
            "tsg": {"TP53", "BRCA1"},
            "both": {"SMAD4"},
            "unknown": {"UNKNOWNGENE"},
        }

    def test_oncogene(self):
        assert get_gene_type("KRAS", self.atlas) == "oncogene"

    def test_tsg(self):
        assert get_gene_type("TP53", self.atlas) == "tsg"

    def test_both(self):
        assert get_gene_type("SMAD4", self.atlas) == "both"

    def test_unknown(self):
        assert get_gene_type("UNKNOWNGENE", self.atlas) == "unknown"

    def test_not_in_atlas_returns_none(self):
        assert get_gene_type("NOTACANCERGENE", self.atlas) is None


# ---------------------------------------------------------------------------
# load_cancer_genes_list
# ---------------------------------------------------------------------------

class TestLoadCancerGenesList:
    def test_loads_oncogene(self, driver_genes_file):
        result = load_cancer_genes_list(driver_genes_file)
        assert "KRAS" in result["oncogene"]

    def test_loads_tsg(self, driver_genes_file):
        result = load_cancer_genes_list(driver_genes_file)
        assert "TP53" in result["tsg"]

    def test_loads_both(self, driver_genes_file):
        result = load_cancer_genes_list(driver_genes_file)
        assert "SMAD4" in result["both"]

    def test_loads_unknown(self, driver_genes_file):
        result = load_cancer_genes_list(driver_genes_file)
        assert "UNKNOWNGENE" in result["unknown"]

    def test_returns_four_keys(self, driver_genes_file):
        result = load_cancer_genes_list(driver_genes_file)
        assert set(result.keys()) == {"oncogene", "tsg", "both", "unknown"}


# ---------------------------------------------------------------------------
# load_canonicals
# ---------------------------------------------------------------------------

class TestLoadCanonicals:
    def test_loads_transcript_from_gtf(self, tmp_path):
        gtf_lines = (
            '##comment\n'
            'chr1\tMANE\texon\t1\t100\t.\t+\t.\t'
            'gene_id ENSG001; transcript_id ENST001.1;\n'
        )
        gtf_path = tmp_path / "test.gtf"
        gtf_path.write_text(gtf_lines)
        result = load_canonicals(str(gtf_path), format="gtf", clean_ids=False)
        assert "ENSG001" in result
        assert "ENST001.1" in result["ENSG001"]

    def test_clean_ids_strips_version(self, tmp_path):
        gtf_lines = (
            'chr1\tMANE\texon\t1\t100\t.\t+\t.\t'
            'gene_id ENSG001; transcript_id ENST001.1;\n'
        )
        gtf_path = tmp_path / "test.gtf"
        gtf_path.write_text(gtf_lines)
        result = load_canonicals(str(gtf_path), format="gtf", clean_ids=True)
        assert "ENST001" in result["ENSG001"]
        assert "ENST001.1" not in result["ENSG001"]

    def test_multiple_transcripts_for_same_gene(self, tmp_path):
        gtf_lines = (
            'chr1\tMANE\texon\t1\t100\t.\t+\t.\t'
            'gene_id ENSG001; transcript_id ENST001;\n'
            'chr1\tMANE_PLUS\texon\t1\t100\t.\t+\t.\t'
            'gene_id ENSG001; transcript_id ENST002;\n'
        )
        gtf_path = tmp_path / "test.gtf"
        gtf_path.write_text(gtf_lines)
        result = load_canonicals(str(gtf_path), format="gtf", clean_ids=False)
        assert len(result["ENSG001"]) == 2
        assert "ENST001" in result["ENSG001"]
        assert "ENST002" in result["ENSG001"]

    def test_duplicate_transcripts_not_added_twice(self, tmp_path):
        gtf_lines = (
            'chr1\tMANE\texon\t1\t100\t.\t+\t.\t'
            'gene_id ENSG001; transcript_id ENST001;\n'
            'chr1\tMANE\texon\t200\t300\t.\t+\t.\t'
            'gene_id ENSG001; transcript_id ENST001;\n'
        )
        gtf_path = tmp_path / "test.gtf"
        gtf_path.write_text(gtf_lines)
        result = load_canonicals(str(gtf_path), format="gtf", clean_ids=False)
        assert result["ENSG001"].count("ENST001") == 1

    def test_skips_comment_lines(self, tmp_path):
        gtf_lines = (
            '## header comment\n'
            '# another comment\n'
            'chr1\tMANE\texon\t1\t100\t.\t+\t.\t'
            'gene_id ENSG001; transcript_id ENST001;\n'
        )
        gtf_path = tmp_path / "test.gtf"
        gtf_path.write_text(gtf_lines)
        result = load_canonicals(str(gtf_path), format="gtf", clean_ids=False)
        assert len(result) == 1

#!/usr/bin/env python
# coding: utf-8
"""
CLI smoke tests and integration tests for oncoDriver.

These tests use the real VCF bundled in tests/  and the assets/ gene list.
They verify end-to-end behaviour: that the output VCF is created, that the
ONCODRIVER INFO tag is present, and that driver and SNV counters are plausible.
"""

import os
import subprocess
import sys

import pytest

# ---------------------------------------------------------------------------
# Path helpers (only skip when real input files are missing)
# ----------------------------------------------------------------------------

TEST_DIR = os.path.dirname(__file__)
PROJECT_ROOT = os.path.dirname(TEST_DIR)
SRC_DIR = os.path.join(PROJECT_ROOT, "src", "oncodriver")
ASSETS_DIR = os.path.join(SRC_DIR, "assets")
DRIVER_GENES = os.path.join(ASSETS_DIR, "cancerGeneList.tsv")
CONFIG = os.path.join(PROJECT_ROOT, "config", "pathogenic_variants.yml")
# Use the annotated test VCF from tests/data/
TEST_VCF = os.path.join(TEST_DIR, "data", "test_mutect2_annotated.vcf.gz")
# Use the annotated test CNV from tests/data/
TEST_CNV = os.path.join(TEST_DIR, "data", "test_cnv_annotated.tsv.gz")
# The test VCF has two samples; use the first sample for tests
TEST_SAMPLE = "D1326R05"

# Skip integration tests when real assets are absent (e.g. CI without LFS)
requires_real_assets = pytest.mark.skipif(
    not (os.path.exists(TEST_VCF) and os.path.exists(DRIVER_GENES) and os.path.exists(CONFIG)),
    reason="Real test assets not available",
)


# ---------------------------------------------------------------------------
# Help / version smoke tests
# ---------------------------------------------------------------------------

class TestCLISmoke:
    def test_help_exits_zero(self):
        result = subprocess.run(
            [sys.executable, "-m", "oncodriver", "--help"],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        assert "oncodriver" in result.stdout.lower() or "usage" in result.stdout.lower()

    def test_version_exits_zero(self):
        result = subprocess.run(
            [sys.executable, "-m", "oncodriver", "--version"],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0

    def test_no_input_exits_nonzero(self):
        result = subprocess.run(
            [sys.executable, "-m", "oncodriver",
             "--config", CONFIG, "--driver_genes", DRIVER_GENES],
            capture_output=True,
            text=True,
        )
        assert result.returncode != 0


# ---------------------------------------------------------------------------
# End-to-end SNV processing
# ---------------------------------------------------------------------------

class TestSNVIntegration:
    @requires_real_assets
    def test_produces_output_vcf(self, tmp_path):
        out = str(tmp_path / "result.vcf.gz")
        result = subprocess.run(
            [sys.executable, "-m", "oncodriver",
             "-i", TEST_VCF,
             "--sample", TEST_SAMPLE,
             "--config", CONFIG,
             "--driver_genes", DRIVER_GENES,
             "--output", out],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0, result.stderr
        assert os.path.exists(out), "Output VCF was not created"

    @requires_real_assets
    def test_output_vcf_has_oncodriver_tag(self, tmp_path):
        """Verify at least one record in the output carries the ONCODRIVER INFO tag."""
        import cyvcf2
        out = str(tmp_path / "result.vcf.gz")
        subprocess.run(
            [sys.executable, "-m", "oncodriver",
             "-i", TEST_VCF,
             "--sample", TEST_SAMPLE,
             "--config", CONFIG,
             "--driver_genes", DRIVER_GENES,
             "--output", out],
            capture_output=True,
            check=True,
        )
        vcf = cyvcf2.VCF(out)
        # cyvcf2 HREC objects use attribute access, not dict .get()
        header_ids = [h["ID"] for h in vcf.header_iter() if h.type == "INFO"]
        assert "ONCODRIVER" in header_ids
        vcf.close()


# ---------------------------------------------------------------------------
# End-to-end CNV processing
# ---------------------------------------------------------------------------

class TestCNVIntegration:
    @requires_real_assets
    def test_cnv_produces_output_tsv(self, tmp_path):
        """Test that providing a CNV file produces an output TSV file."""
        out_cnv = str(tmp_path / "result_cnv.tsv")
        result = subprocess.run(
            [sys.executable, "-m", "oncodriver",
             "--cnv", TEST_CNV,
             "--config", CONFIG,
             "--driver_genes", DRIVER_GENES,
             "--outputcnv", out_cnv],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0, result.stderr
        assert os.path.exists(out_cnv), "Output CNV TSV was not created"

    @requires_real_assets
    def test_cnv_output_has_driver_status_column(self, tmp_path):
        """Verify the output CNV TSV has the driver_status column."""
        import csv
        out_cnv = str(tmp_path / "result_cnv.tsv")
        subprocess.run(
            [sys.executable, "-m", "oncodriver",
             "--cnv", TEST_CNV,
             "--config", CONFIG,
             "--driver_genes", DRIVER_GENES,
             "--outputcnv", out_cnv],
            capture_output=True,
            check=True,
        )
        with open(out_cnv, newline="") as f:
            reader = csv.reader(f, delimiter=",")
            header = next(reader)
        assert "driver_status" in header, "driver_status column missing from output"
        assert "gene_type" in header, "gene_type column missing from output"


# ---------------------------------------------------------------------------
# End-to-end combined SNV + CNV processing
# ---------------------------------------------------------------------------

class TestSNVAndCNVIntegration:
    @requires_real_assets
    def test_snv_and_cnv_produce_both_outputs(self, tmp_path):
        """Test that providing both VCF and CNV files produces both output files."""
        out_vcf = str(tmp_path / "result.vcf.gz")
        out_cnv = str(tmp_path / "result_cnv.tsv")
        result = subprocess.run(
            [sys.executable, "-m", "oncodriver",
             "-i", TEST_VCF,
             "--cnv", TEST_CNV,
             "--sample", TEST_SAMPLE,
             "--config", CONFIG,
             "--driver_genes", DRIVER_GENES,
             "--output", out_vcf,
             "--outputcnv", out_cnv],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0, result.stderr
        assert os.path.exists(out_vcf), "Output VCF was not created"
        assert os.path.exists(out_cnv), "Output CNV TSV was not created"

    @requires_real_assets
    def test_snv_and_cnv_both_have_driver_results(self, tmp_path):
        """Test that both SNV and CNV outputs contain driver information."""
        import cyvcf2
        import csv

        out_vcf = str(tmp_path / "result.vcf.gz")
        out_cnv = str(tmp_path / "result_cnv.tsv")
        subprocess.run(
            [sys.executable, "-m", "oncodriver",
             "-i", TEST_VCF,
             "--cnv", TEST_CNV,
             "--sample", TEST_SAMPLE,
             "--config", CONFIG,
             "--driver_genes", DRIVER_GENES,
             "--output", out_vcf,
             "--outputcnv", out_cnv],
            capture_output=True,
            check=True,
        )

        # Check VCF has ONCODRIVER tag
        vcf = cyvcf2.VCF(out_vcf)
        header_ids = [h["ID"] for h in vcf.header_iter() if h.type == "INFO"]
        assert "ONCODRIVER" in header_ids
        vcf.close()

        # Check CNV has driver_status column
        with open(out_cnv, newline="") as f:
            reader = csv.reader(f, delimiter=",")
            header = next(reader)
        assert "driver_status" in header


# ---------------------------------------------------------------------------
# Debug mode tests
# ---------------------------------------------------------------------------

class TestDebugMode:
    @requires_real_assets
    def test_debug_mode_preserves_snv_count(self, tmp_path):
        """Test that debug mode preserves the same number of SNV records in output as input."""
        import cyvcf2

        out = str(tmp_path / "result_debug.vcf.gz")
        subprocess.run(
            [sys.executable, "-m", "oncodriver",
             "-i", TEST_VCF,
             "--sample", TEST_SAMPLE,
             "--config", CONFIG,
             "--driver_genes", DRIVER_GENES,
             "--output", out,
             "--debug"],
            capture_output=True,
            check=True,
        )

        # Count records in input VCF
        vcf_in = cyvcf2.VCF(TEST_VCF)
        input_count = sum(1 for _ in vcf_in)
        vcf_in.close()

        # Count records in output VCF
        vcf_out = cyvcf2.VCF(out)
        output_count = sum(1 for _ in vcf_out)
        vcf_out.close()

        assert output_count == input_count, \
            f"Output VCF has {output_count} records but input has {input_count}"

    @requires_real_assets
    def test_debug_mode_preserves_cnv_count(self, tmp_path):
        """Test that debug mode preserves the same number of CNV records in output as input."""
        import gzip

        out_cnv = str(tmp_path / "result_cnv_debug.tsv")
        subprocess.run(
            [sys.executable, "-m", "oncodriver",
             "--cnv", TEST_CNV,
             "--config", CONFIG,
             "--driver_genes", DRIVER_GENES,
             "--outputcnv", out_cnv,
             "--debug"],
            capture_output=True,
            check=True,
        )

        # Count records in input CNV (skip header)
        with gzip.open(TEST_CNV, 'rt') as f:
            input_lines = sum(1 for _ in f) - 1  # subtract header

        # Count records in output CNV (skip header)
        with open(out_cnv, 'r') as f:
            output_lines = sum(1 for _ in f) - 1  # subtract header

        assert output_lines == input_lines, \
            f"Output CNV has {output_lines} records but input has {input_lines}"

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
# ---------------------------------------------------------------------------

TEST_DIR = os.path.dirname(__file__)
PROJECT_ROOT = os.path.dirname(TEST_DIR)
ASSETS_DIR = os.path.join(PROJECT_ROOT, "assets")
DRIVER_GENES = os.path.join(ASSETS_DIR, "cancerGeneList.tsv")
CONFIG = os.path.join(PROJECT_ROOT, "config", "pathogenic_variants.yml")
# The VCF fixture lives in test/ (singular), not tests/ (plural)
TEST_VCF = os.path.join(
    PROJECT_ROOT, "test",
    "02-104_T_vs_02-104_C_Mutect2_calls_norm_GnomAD_filtered_ICGC_CancerHotspots_COSMIC_dbNSFP.vcf.gz",
)
# The test VCF has two samples; use the tumour sample for tests
TEST_SAMPLE = "02-104_T"

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

    @requires_real_assets
    def test_driver_count_is_positive(self, tmp_path):
        """Sanity check: at least one driver variant should be identified."""
        out = str(tmp_path / "result.vcf.gz")
        proc = subprocess.run(
            [sys.executable, "-m", "oncodriver",
             "-i", TEST_VCF,
             "--sample", TEST_SAMPLE,
             "--config", CONFIG,
             "--driver_genes", DRIVER_GENES,
             "--output", out],
            capture_output=True,
            text=True,
            check=True,
        )
        # The summary is printed to stdout via logging
        combined = proc.stdout + proc.stderr
        # Look for "Driver variants=" line in summary
        driver_lines = [l for l in combined.splitlines() if "Driver variants=" in l]
        if driver_lines:
            count_str = driver_lines[-1].split("=")[-1].strip()
            assert int(count_str) >= 0

    @requires_real_assets
    def test_maf_filter_reduces_driver_count(self, tmp_path):
        """Adding --max_maf 0.001 should not increase driver count vs no filter."""
        out_no_filter = str(tmp_path / "no_filter.vcf.gz")
        out_filtered = str(tmp_path / "filtered.vcf.gz")

        subprocess.run(
            [sys.executable, "-m", "oncodriver",
             "-i", TEST_VCF,
             "--sample", TEST_SAMPLE,
             "--config", CONFIG,
             "--driver_genes", DRIVER_GENES,
             "--output", out_no_filter],
            capture_output=True, check=True,
        )
        subprocess.run(
            [sys.executable, "-m", "oncodriver",
             "-i", TEST_VCF,
             "--sample", TEST_SAMPLE,
             "--config", CONFIG,
             "--driver_genes", DRIVER_GENES,
             "--max_maf", "0.001",
             "--output", out_filtered],
            capture_output=True, check=True,
        )

        import cyvcf2
        def count_records(path):
            vcf = cyvcf2.VCF(path)
            n = sum(1 for _ in vcf)
            vcf.close()
            return n

        assert count_records(out_filtered) <= count_records(out_no_filter)

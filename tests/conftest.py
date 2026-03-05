#!/usr/bin/env python
# coding: utf-8
"""
Shared pytest fixtures for oncoDriver test suite.
"""

import csv
import os
import textwrap
import pytest


# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

TEST_DIR = os.path.dirname(__file__)
ASSETS_DIR = os.path.join(os.path.dirname(TEST_DIR), "assets")
TEST_VCF = os.path.join(TEST_DIR, "02-104_T_vs_02-104_C_Mutect2_calls_norm_GnomAD_filtered_ICGC_CancerHotspots_COSMIC_dbNSFP.vcf.gz")
DRIVER_GENES_FILE = os.path.join(ASSETS_DIR, "cancerGeneList.tsv")
CONFIG_FILE = os.path.join(os.path.dirname(TEST_DIR), "config", "pathogenic_variants.yml")


# ---------------------------------------------------------------------------
# Small in-memory fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def cancer_atlas():
    """Minimal cancer atlas with one gene in each category."""
    return {
        "oncogene": {"KRAS", "EGFR"},
        "tsg": {"TP53", "BRCA1"},
        "both": {"SMAD4"},
        "unknown": {"UNKNOWNGENE"},
    }


@pytest.fixture
def driver_genes_file(tmp_path):
    """Write a minimal oncoKB-format TSV driver genes file."""
    tsv_path = tmp_path / "cancerGeneList.tsv"
    rows = [
        ["Hugo Symbol", "Is Oncogene", "Is Tumor Suppressor Gene"],
        ["KRAS", "Yes", "No"],
        ["TP53", "No", "Yes"],
        ["SMAD4", "Yes", "Yes"],
        ["UNKNOWNGENE", "No", "No"],
    ]
    with open(tsv_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerows(rows)
    return str(tsv_path)


@pytest.fixture
def minimal_config_file(tmp_path):
    """Write a minimal valid YAML config file."""
    config_path = tmp_path / "config.yml"
    config_path.write_text(textwrap.dedent("""\
        vcf:
          tag: 'ANN'
          sep: '&'
          vaf: 'AF'
          depth: 'DP'
          alt_depth: 'AD'
          annot_info: 1
          gene_id: 3
          transcript_id: 6
          cancer_db:
            - 'CancerHotspots'
          polym_db:
            gnomAD:
              - AF

        cnv:
          gene_id: 10
          class: 8
          start: 1
          end: 2

        select:
          - var_type: 'snv'
            gene_type: 'oncogene|tsg|both|unknown'
            var_classes:
              - missense_variant
            is_hotspot: true
          - var_type: 'snv'
            gene_type: 'tsg|both|unknown'
            var_classes:
              - stop_gained
            is_hotspot: false
          - var_type: 'cnv'
            gene_type: 'oncogene|both'
            var_classes:
              - AMP
            maxsize: 10000000
          - var_type: 'cnv'
            gene_type: 'tsg|both'
            var_classes:
              - DEL
    """))
    return str(config_path)


@pytest.fixture
def minimal_cnv_file(tmp_path):
    """Write a minimal CNV TSV file (11 columns, gene in column 10)."""
    cnv_path = tmp_path / "segments.tsv"
    header = ["chrom", "loc.start", "loc.end", "ID", "CNt", "Geno",
              "logratio", "ploidy", "call", "LOH", "gene"]
    rows = [
        ["chr1", "1000", "2000", "S1", "8", "AA", "2.0", "2", "AMP", "0", "KRAS"],
        ["chr17", "5000", "6000", "S1", "0", "Del", "-2.0", "2", "DEL", "1", "TP53"],
        ["chr7", "1000", "20000000", "S1", "10", "AA", "3.0", "2", "AMP", "0", "EGFR"],
        ["chr1", "9000", "9100", "S1", "4", "Aa", "0.5", "2", "GAIN", "0", "NOTCANCER"],
    ]
    with open(cnv_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(header)
        writer.writerows(rows)
    return str(cnv_path)

#!/usr/bin/env python
# coding: utf-8
"""
Unit and integration tests for oncodriver.cnv_process module.
"""

import csv
import os
import pytest
from oncodriver.cnv_process import is_driver_cnv, process_cnv


# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

TEST_DIR = os.path.dirname(__file__)
TEST_CNV = os.path.join(TEST_DIR, "data", "test_cnv_annotated.tsv.gz")


# ---------------------------------------------------------------------------
# is_driver_cnv
# ---------------------------------------------------------------------------

class TestIsDriverCnv:
    """Tests for is_driver_cnv() – critical due to the loop-exit bug (fix 1.1)."""

    @pytest.fixture
    def amp_rule_with_maxsize(self):
        return [{"var_classes": ["AMP"], "maxsize": 10_000_000}]

    @pytest.fixture
    def del_rule(self):
        return [{"var_classes": ["DEL"]}]

    @pytest.fixture
    def two_rules(self):
        """Two rules to verify all rules are evaluated (regression for fix 1.1)."""
        return [
            {"var_classes": ["AMP"], "maxsize": 10_000_000},
            {"var_classes": ["DEL"]},
        ]

    # AMP tests
    def test_amp_below_maxsize_is_driver(self, amp_rule_with_maxsize):
        assert is_driver_cnv("AMP", 5_000_000, amp_rule_with_maxsize) is True

    def test_amp_above_maxsize_is_not_driver(self, amp_rule_with_maxsize):
        assert is_driver_cnv("AMP", 15_000_000, amp_rule_with_maxsize) is False

    def test_amp_equal_maxsize_is_driver(self, amp_rule_with_maxsize):
        # length == maxsize: condition is length > maxsize, so equal is still a driver
        assert is_driver_cnv("AMP", 10_000_000, amp_rule_with_maxsize) is True

    # DEL tests
    def test_del_is_driver(self, del_rule):
        assert is_driver_cnv("DEL", 1_000_000, del_rule) is True

    def test_gain_not_driver_against_del_rule(self, del_rule):
        assert is_driver_cnv("GAIN", 1_000_000, del_rule) is False

    # Regression: second rule must be evaluated (fix 1.1)
    def test_second_rule_evaluated_when_first_doesnt_match(self, two_rules):
        # First rule: AMP → no match for DEL
        # Second rule: DEL → match
        assert is_driver_cnv("DEL", 500_000, two_rules) is True

    def test_unknown_class_returns_false(self, amp_rule_with_maxsize):
        assert is_driver_cnv("UNKNOWN", 1_000_000, amp_rule_with_maxsize) is False

    def test_amp_no_maxsize_rule_is_driver(self):
        rules = [{"var_classes": ["AMP"]}]
        assert is_driver_cnv("AMP", 999_999_999, rules) is True


# ---------------------------------------------------------------------------
# process_cnv – integration smoke tests
# ---------------------------------------------------------------------------

class TestProcessCnv:
    @pytest.fixture
    def cnv_conf(self):
        return {"gene_id": 10, "class": 8, "start": 1, "end": 2}

    @pytest.fixture
    def select_algo(self):
        return {
            ("oncogene", "cnv"): [{"var_classes": ["AMP"], "maxsize": 10_000_000}],
            ("tsg", "cnv"): [{"var_classes": ["DEL"]}],
        }

    @pytest.fixture
    def driver_genes(self):
        return {
            "oncogene": {"KRAS", "EGFR"},
            "tsg": {"TP53"},
            "both": set(),
            "unknown": set(),
        }

    def test_amp_in_oncogene_is_driver(self, tmp_path, minimal_cnv_file, cnv_conf, select_algo, driver_genes):
        out = str(tmp_path / "out.tsv")
        driver_count, cnv_count = process_cnv(
            minimal_cnv_file, out, debug=False,
            cnv_conf=cnv_conf, select_algo=select_algo,
            driver_genes=driver_genes, can_ids=None,
        )
        assert cnv_count == 4  # 4 rows in fixture
        # KRAS (oncogene) AMP → driver
        # TP53 (tsg) DEL → driver
        assert driver_count == 2

    def test_output_file_created(self, tmp_path, minimal_cnv_file, cnv_conf, select_algo, driver_genes):
        out = str(tmp_path / "out.tsv")
        process_cnv(minimal_cnv_file, out, debug=False,
                    cnv_conf=cnv_conf, select_algo=select_algo,
                    driver_genes=driver_genes, can_ids=None)
        assert os.path.exists(out)

    def test_output_has_correct_header(self, tmp_path, minimal_cnv_file, cnv_conf, select_algo, driver_genes):
        out = str(tmp_path / "out.tsv")
        process_cnv(minimal_cnv_file, out, debug=False,
                    cnv_conf=cnv_conf, select_algo=select_algo,
                    driver_genes=driver_genes, can_ids=None)
        with open(out, newline="") as f:
            # csv.writer uses comma by default; read with matching delimiter
            reader = csv.reader(f)
            header = next(reader)
        assert "driver_status" in header
        assert "gene_type" in header

    def test_debug_mode_writes_non_drivers_but_does_not_count_them(
        self, tmp_path, minimal_cnv_file, cnv_conf, select_algo, driver_genes
    ):
        """Regression test for fix 1.2: driver_counter must not include non-driver rows."""
        out = str(tmp_path / "out_debug.tsv")
        driver_count, cnv_count = process_cnv(
            minimal_cnv_file, out, debug=True,
            cnv_conf=cnv_conf, select_algo=select_algo,
            driver_genes=driver_genes, can_ids=None,
        )
        # Even in debug mode, driver_count should only reflect true drivers
        assert driver_count == 2

    def test_cnv_returns_zero_when_no_file(self, cnv_conf, select_algo, driver_genes):
        driver_count, cnv_count = process_cnv(
            None, None, debug=False,
            cnv_conf=cnv_conf, select_algo=select_algo,
            driver_genes=driver_genes, can_ids=None,
        )
        assert driver_count == 0
        assert cnv_count == 0

    def test_gene_id_rule_overrides_gene_type(self, tmp_path, cnv_conf, driver_genes):
        """Regression for fix 1.3: gene-id rule should work even if gene_type is None."""
        # SPECIALGENE is not in driver_genes, but has a gene_id rule
        select_with_gene_id = {
            ("SPECIALGENE", "cnv"): [{"var_classes": ["AMP"]}],
        }
        cnv_path = tmp_path / "special.tsv"
        header = ["chrom", "loc.start", "loc.end", "ID", "CNt", "Geno",
                  "logratio", "ploidy", "call", "LOH", "gene"]
        with open(cnv_path, "w", newline="") as f:
            writer = csv.writer(f, delimiter="\t")
            writer.writerow(header)
            writer.writerow(["chr1", "100", "200", "S1", "8", "AA",
                              "2.0", "2", "AMP", "0", "SPECIALGENE"])
        out = str(tmp_path / "special_out.tsv")
        driver_count, cnv_count = process_cnv(
            str(cnv_path), out, debug=False,
            cnv_conf=cnv_conf, select_algo=select_with_gene_id,
            driver_genes=driver_genes, can_ids=None,
        )
        assert driver_count == 1


# ---------------------------------------------------------------------------
# Integration tests with annotated CNV file
# ---------------------------------------------------------------------------

class TestProcessCnvAnnotated:
    """Integration tests using the annotated CNV test file."""

    @pytest.fixture
    def cnv_conf(self):
        return {"gene_id": 10, "class": 8, "start": 1, "end": 2}

    @pytest.fixture
    def select_algo(self):
        return {
            ("oncogene", "cnv"): [{"var_classes": ["AMP"], "maxsize": 10_000_000}],
            ("tsg", "cnv"): [{"var_classes": ["DEL"]}],
        }

    @pytest.fixture
    def driver_genes(self):
        return {
            "oncogene": {"KRAS", "EGFR"},
            "tsg": {"TP53"},
            "both": set(),
            "unknown": set(),
        }

    def test_annotated_cnv_file_exists(self):
        """Verify the annotated CNV test file exists."""
        assert os.path.exists(TEST_CNV)

    def test_process_annotated_cnv(self, tmp_path, cnv_conf, select_algo, driver_genes):
        """Test processing the annotated CNV file."""
        out = str(tmp_path / "annotated_out.tsv")
        driver_count, cnv_count = process_cnv(
            TEST_CNV, out, debug=False,
            cnv_conf=cnv_conf, select_algo=select_algo,
            driver_genes=driver_genes, can_ids=None,
        )
        # Should process all rows from the annotated file
        assert cnv_count > 0
        # Check output file was created
        assert os.path.exists(out)

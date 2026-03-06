#!/usr/bin/env python
# coding: utf-8
"""
Unit tests for oncodriver.utils module.
"""

import textwrap
import pytest
from oncodriver.utils import load_configuration, print_config


# ---------------------------------------------------------------------------
# load_configuration
# ---------------------------------------------------------------------------

class TestLoadConfiguration:
    def test_loads_valid_config(self, minimal_config_file):
        cfg = load_configuration(minimal_config_file)
        assert "vcf" in cfg
        assert "cnv" in cfg
        assert "select" in cfg

    def test_select_contains_gene_type_keys(self, minimal_config_file):
        cfg = load_configuration(minimal_config_file)
        # oncogene + snv rule should exist after parsing
        assert ("oncogene", "snv") in cfg["select"]

    def test_select_splits_pipe_separated_gene_types(self, minimal_config_file):
        cfg = load_configuration(minimal_config_file)
        # 'tsg|both|unknown' → three separate keys
        for gene_type in ("tsg", "both", "unknown"):
            assert (gene_type, "snv") in cfg["select"]

    def test_missing_vcf_section_raises(self, tmp_path):
        cfg_path = tmp_path / "bad.yml"
        cfg_path.write_text(textwrap.dedent("""\
            cnv:
              gene_id: 10
            select: []
        """))
        with pytest.raises(ValueError, match="vcf"):
            load_configuration(str(cfg_path))

    def test_missing_cnv_section_raises(self, tmp_path):
        cfg_path = tmp_path / "bad.yml"
        cfg_path.write_text(textwrap.dedent("""\
            vcf:
              tag: ANN
            select: []
        """))
        with pytest.raises(ValueError, match="cnv"):
            load_configuration(str(cfg_path))

    def test_missing_select_section_raises(self, tmp_path):
        cfg_path = tmp_path / "bad.yml"
        cfg_path.write_text(textwrap.dedent("""\
            vcf:
              tag: ANN
            cnv:
              gene_id: 10
        """))
        with pytest.raises(ValueError, match="select"):
            load_configuration(str(cfg_path))

    def test_malformed_yaml_raises_value_error(self, tmp_path):
        cfg_path = tmp_path / "bad.yml"
        cfg_path.write_text("key: [unclosed bracket\n")
        with pytest.raises(ValueError):
            load_configuration(str(cfg_path))

    def test_gene_id_rule_creates_correct_key(self, tmp_path):
        cfg_path = tmp_path / "gene_id.yml"
        cfg_path.write_text(textwrap.dedent("""\
            vcf:
              tag: ANN
            cnv:
              gene_id: 10
            select:
              - var_type: snv
                gene_id: 'TERT'
                var_classes:
                  - upstream_gene_variant
        """))
        cfg = load_configuration(str(cfg_path))
        assert ("TERT", "snv") in cfg["select"]

    def test_pipe_separated_gene_ids_create_multiple_keys(self, tmp_path):
        cfg_path = tmp_path / "pipe_ids.yml"
        cfg_path.write_text(textwrap.dedent("""\
            vcf:
              tag: ANN
            cnv:
              gene_id: 10
            select:
              - var_type: cnv
                gene_id: 'SOX2|TP63'
                var_classes:
                  - AMP
        """))
        cfg = load_configuration(str(cfg_path))
        assert ("SOX2", "cnv") in cfg["select"]
        assert ("TP63", "cnv") in cfg["select"]


# ---------------------------------------------------------------------------
# print_config
# ---------------------------------------------------------------------------

class TestPrintConfig:
    def test_runs_without_error_flat_dict(self, caplog):
        import logging
        with caplog.at_level(logging.INFO):
            print_config({"version": "1.0", "mode": "strict"})
        assert "VERSION" in caplog.text

    def test_runs_without_error_nested_dict(self, caplog):
        import logging
        cfg = {"vcf": {"tag": "ANN", "vaf": "AF"}}
        with caplog.at_level(logging.INFO):
            print_config(cfg)
        assert "VCF" in caplog.text
        assert "tag" in caplog.text

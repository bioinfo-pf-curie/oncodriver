#!/usr/bin/env python
# coding: utf-8
"""
Unit tests for oncodriver.snv_process module.
"""

import numpy as np
import pytest
from oncodriver.snv_process import (
    get_INFO,
    subset_INFO,
    is_hotspot,
    is_polym,
    is_driver,
)


# ---------------------------------------------------------------------------
# get_INFO
# ---------------------------------------------------------------------------

class TestGetINFO:
    def test_parses_single_annotation(self):
        ann = "T|missense_variant|MODERATE|KRAS|ENSG00001|transcript|ENST00001|protein_coding"
        result = get_INFO(ann)
        assert len(result) == 1
        assert result[0][1] == "missense_variant"

    def test_parses_multiple_annotations(self):
        ann = "T|missense_variant|MODERATE|KRAS|ENS1|t|ENST1|pc,T|stop_gained|HIGH|TP53|ENS2|t|ENST2|pc"
        result = get_INFO(ann)
        assert len(result) == 2
        assert result[0][1] == "missense_variant"
        assert result[1][1] == "stop_gained"

    def test_returns_empty_list_for_none(self):
        assert get_INFO(None) == []

    def test_returns_empty_list_for_empty_string(self):
        assert get_INFO("") == []

    def test_returns_list_of_dicts(self):
        ann = "A|B|C"
        result = get_INFO(ann)
        assert isinstance(result, list)
        assert isinstance(result[0], dict)


# ---------------------------------------------------------------------------
# subset_INFO
# ---------------------------------------------------------------------------

class TestSubsetINFO:
    def test_subsets_dict(self):
        annot = {"AF": 0.05, "DP": 100, "OTHER": "x"}
        result = subset_INFO(annot, keys=["AF", "DP"])
        assert result == {"AF": 0.05, "DP": 100}

    def test_missing_keys_ignored(self):
        annot = {"AF": 0.05}
        result = subset_INFO(annot, keys=["AF", "NOTPRESENT"])
        assert result == {"AF": 0.05}

    def test_subsets_list_of_dicts(self):
        annot = [{"AF": 0.05, "X": 1}, {"AF": 0.1, "Y": 2}]
        result = subset_INFO(annot, keys=["AF"])
        assert len(result) == 2
        assert result[0] == {"AF": 0.05}

    def test_empty_keys_returns_empty(self):
        annot = {"AF": 0.05, "DP": 100}
        result = subset_INFO(annot, keys=[])
        assert result == {}


# ---------------------------------------------------------------------------
# is_hotspot
# ---------------------------------------------------------------------------

class TestIsHotspot:
    def test_flag_present_returns_true(self):
        infos = {"CancerHotspots": "1", "OTHER": "X"}
        assert is_hotspot(None, infos=infos, flags=["CancerHotspots"]) is True

    def test_flag_dot_returns_false(self):
        infos = {"CancerHotspots": "."}
        assert is_hotspot(None, infos=infos, flags=["CancerHotspots"]) is False

    def test_flag_none_returns_false(self):
        infos = {"CancerHotspots": None}
        assert is_hotspot(None, infos=infos, flags=["CancerHotspots"]) is False

    def test_flag_missing_returns_false(self):
        infos = {"OTHER": "1"}
        assert is_hotspot(None, infos=infos, flags=["CancerHotspots"]) is False

    def test_empty_flags_returns_false(self):
        infos = {"CancerHotspots": "1"}
        assert is_hotspot(None, infos=infos, flags=[]) is False


# ---------------------------------------------------------------------------
# is_polym
# ---------------------------------------------------------------------------

class TestIsPolym:
    @pytest.fixture
    def flags(self):
        return {"gnomAD": ["AF"]}

    def test_maf_above_threshold_returns_true(self, flags):
        infos = {"AF": 0.05}
        assert is_polym(None, infos=infos, flags=flags, val=0.01) is True

    def test_maf_equal_threshold_returns_true(self, flags):
        infos = {"AF": 0.01}
        assert is_polym(None, infos=infos, flags=flags, val=0.01) is True

    def test_maf_below_threshold_returns_false(self, flags):
        infos = {"AF": 0.001}
        assert is_polym(None, infos=infos, flags=flags, val=0.01) is False

    def test_maf_dot_returns_false(self, flags):
        infos = {"AF": "."}
        assert is_polym(None, infos=infos, flags=flags, val=0.01) is False

    def test_maf_none_returns_false(self, flags):
        infos = {"AF": None}
        assert is_polym(None, infos=infos, flags=flags, val=0.01) is False

    def test_maf_na_returns_false(self, flags):
        infos = {"AF": "NA"}
        assert is_polym(None, infos=infos, flags=flags, val=0.01) is False

    def test_tuple_value_above_threshold(self, flags):
        infos = {"AF": (0.001, 0.05)}
        assert is_polym(None, infos=infos, flags=flags, val=0.01) is True

    def test_tuple_value_all_below_threshold(self, flags):
        infos = {"AF": (0.001, 0.0001)}
        assert is_polym(None, infos=infos, flags=flags, val=0.01) is False


# ---------------------------------------------------------------------------
# is_driver (SNV)
# ---------------------------------------------------------------------------

class TestIsDriver:
    """Tests for the is_driver() function that applies the selection rules."""

    @pytest.fixture
    def vcf_conf(self):
        return {
            "annot_info": 1,
            "sep": "&",
            "cancer_db": ["CancerHotspots"],
        }

    @pytest.fixture
    def conf_select_hotspot_required(self):
        return [{"var_classes": ["missense_variant"], "is_hotspot": True}]

    @pytest.fixture
    def conf_select_no_hotspot_required(self):
        return [{"var_classes": ["stop_gained"], "is_hotspot": False}]

    def test_driver_when_hotspot_and_class_match(self, vcf_conf, conf_select_hotspot_required):
        variant = {1: "missense_variant"}
        db_info = {"CancerHotspots": "1"}
        result, decision = is_driver(variant, db_info, conf_select_hotspot_required, vcf_conf)
        assert result is True
        assert decision == "missense_variant"

    def test_not_driver_when_hotspot_missing(self, vcf_conf, conf_select_hotspot_required):
        variant = {1: "missense_variant"}
        db_info = {"CancerHotspots": "."}
        result, decision = is_driver(variant, db_info, conf_select_hotspot_required, vcf_conf)
        assert result is False
        assert decision is None

    def test_driver_when_hotspot_not_required(self, vcf_conf, conf_select_no_hotspot_required):
        variant = {1: "stop_gained"}
        db_info = {}
        result, decision = is_driver(variant, db_info, conf_select_no_hotspot_required, vcf_conf)
        assert result is True
        assert decision == "stop_gained"

    def test_not_driver_when_var_class_not_in_rules(self, vcf_conf, conf_select_no_hotspot_required):
        variant = {1: "synonymous_variant"}
        db_info = {}
        result, decision = is_driver(variant, db_info, conf_select_no_hotspot_required, vcf_conf)
        assert result is False
        assert decision is None

    def test_sep_splits_annotation(self, vcf_conf, conf_select_no_hotspot_required):
        """Test that pipe/amp-separated annotations are split correctly."""
        variant = {1: "synonymous_variant&stop_gained"}
        db_info = {}
        result, decision = is_driver(variant, db_info, conf_select_no_hotspot_required, vcf_conf)
        assert result is True

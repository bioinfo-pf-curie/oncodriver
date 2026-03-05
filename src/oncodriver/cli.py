#!/usr/bin/env python
# coding: utf-8
#
#  This file is part of oncoDriver software.
#
#  Copyright (c) 2026 - Institut Curie
#
#  File author(s):
#      Nicolas Servant <nicolas.servant@curie.fr>
#
#  Distributed under the terms of the CeCILL-B license.
#  The full license is in the LICENSE file, distributed with this software.
#
##############################################################################

"""
This script is designed to apply a decision tree on a VCF file
and to select driver mutations from a VCF file.

python -m oncodriver -i ${VCF} \
  --config ${DECISION_TREE} \
  --use_canonical \
  --canonical_reference ${TABLE}  > oncoDriver_results.log
"""

from datetime import date
import argparse
import sys
import os
import logging

from .utils import load_configuration, print_config

logger = logging.getLogger(__name__)
from .annot import load_cancer_genes_list, load_canonicals
from .snv_process import process_vcf
from .cnv_process import process_cnv
from . import __version__


def args_parse() -> argparse.Namespace:
    """
    Parse command-line arguments.

    Returns:
        Parsed arguments namespace.
    """
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", "--vcf", help="Input file (.vcf, .vcf.gz, .bcf)")
    parser.add_argument("--cnv", help="Input file for CNV (segments.transformed_annot_oncokb.txt)", default=None)

    # Configs
    parser.add_argument("--config", help="Decision tree config file", type=str, required=True)
    parser.add_argument("--sample", help="Specify the sample ID to focus on", type=str, default=None)

    # Technical filters
    parser.add_argument("--min_vaf", help="Filter variants with Allelic Ratio <= vaf", type=float, default=None)
    parser.add_argument("--max_maf", help="Filter variants with MAF > maf", type=float, default=None)
    parser.add_argument("--min_depth", help="Filter variants with depth < min_depth", type=int, default=None)
    parser.add_argument("--min_alt_depth", help="Filter variants with alternative allele depth <= min_alt_depth", type=int, default=None)

    # Canonical transcript
    parser.add_argument("--use_canonical", help="Use canonical transcript only", action="store_true")
    parser.add_argument("--canonical_ids", help="File with canonical IDs information", type=str, default=None)
    
    # Driver genes
    parser.add_argument("--driver_genes", help="Cancer gene list from oncoKB", type=str, required=True)
    
    # Others
    parser.add_argument("--output", help="Output file name", type=str)
    parser.add_argument("--outputcnv", help="CNV output file name", type=str)
    parser.add_argument("--verbose", help="Active verbose mode", action="store_true")
    parser.add_argument("--debug", help="Export original VCF with ONCODRIVER tag", action="store_true")
    parser.add_argument("--strict", help="Raise an error instead of warnings", action="store_true")
    parser.add_argument("--version", help="Version number", action='version', version="%(prog)s (" + __version__ + ")")

    args = parser.parse_args()
    return (args)


def parse_and_validate_args() -> argparse.Namespace:
    """
    Parse and validate command-line arguments.
    """
    args = args_parse()
    if args.vcf is None and args.cnv is None:
        logger.error("Please provide at least a vcf ('-i') file or a cnv ('--cnv') file")
        sys.exit(-1)
    if args.use_canonical and args.canonical_ids is None:
        logger.error("Please provide a list of canonical ids with '--canonical_ids'")
        sys.exit(-1)
    return args


def output_summary(vcf_file: str, sample: str, cnv_file: str, config_file: str, driver_genes_file: str, use_canonical: bool, canonical_ids_file: str, conf: dict, snv_counter: int, cnv_counter: int, driver_counter: int) -> None:
    """
    Output the summary of the analysis.
    """
    summary_lines = [
        f"## oncoDriver version={__version__}",
        f"## When={date.today()}",
        f"## SNVs={vcf_file}",
        f"## Sample={sample}",
        f"## CNVs={cnv_file}",
        f"## Config={config_file}",
        f"## Driver genes={driver_genes_file}",
        f"## Use canonical={use_canonical}",
        f"## Canonical genes={canonical_ids_file}",
    ]
    for line in summary_lines:
        logger.info(line)
    print_config(conf)
    logger.info(f"## Total SNVs variants={snv_counter}")
    logger.info(f"## Total CNVs variants={cnv_counter}")
    logger.info(f"## Driver variants={driver_counter}")


def main() -> None:
    """
    Main entry point for the oncoDriver CLI.
    """
    # Configure logging before argument validation so that error messages during
    # validation are properly formatted and visible.
    _level = logging.DEBUG if '--debug' in sys.argv else logging.INFO
    logging.basicConfig(level=_level, format='%(message)s')

    args = parse_and_validate_args()

    conf = load_configuration(args.config)
    select_algo = conf['select']
    vcf_conf = conf['vcf']
    cnv_conf = conf['cnv']

    driver_genes = load_cancer_genes_list(args.driver_genes)
    can_ids = None
    if args.use_canonical:
        can_ids = load_canonicals(args.canonical_ids, clean_ids=True)
    
    driver_counter, snv_counter = process_vcf(args.vcf, args.sample, 
                                              args.output, args.debug, 
                                              args.min_vaf, args.min_depth, 
                                              args.min_alt_depth, args.max_maf, 
                                              args.use_canonical, args.strict, 
                                              vcf_conf, select_algo, driver_genes, can_ids)
    
    driver_counter, cnv_counter = process_cnv(args.cnv, args.outputcnv, 
                                              args.debug, cnv_conf, 
                                              select_algo, driver_genes, 
                                              can_ids, driver_counter)

    output_summary(args.vcf, args.sample, 
                   args.cnv, args.config, 
                   args.driver_genes, args.use_canonical, 
                   args.canonical_ids, conf, 
                   snv_counter, cnv_counter, driver_counter)


if __name__ == "__main__":
    main()

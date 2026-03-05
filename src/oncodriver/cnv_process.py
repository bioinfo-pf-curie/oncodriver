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

import logging
import re
import csv
import gzip
import os
import sys

from .annot import get_gene_type

logger = logging.getLogger(__name__)

DEFAULT_CNV_SUFFIX = '_oncoDriver.tsv'


def is_driver_cnv(var_class: str, length: int, conf_select: dict) -> bool:
    """
    Determine if a copy number variant is driver.
    Args:
        var_class: CNV class (e.g., 'AMP').
        length: Length of the CNV.
        conf_select: Selection rules.
    Returns:
        True if driver, False otherwise.
    """
    for r in conf_select:
        if var_class in r['var_classes']:
            if (var_class == "AMP" and r.get('maxsize') is not None) and (length > r['maxsize']):
                return False
            return True
    return False
    
def process_cnv(cnv_file: str, outputcnv: str, debug: bool, cnv_conf: dict, select_algo: dict, driver_genes: dict,
                can_ids: dict, initial_driver_counter: int) -> tuple[int, int]:
    """
    Process CNV file and return updated driver and CNV counters.
    Args:
        cnv_file: Path to the CNV file.
        outputcnv: Output CNV file name.
        debug: Debug mode flag.
        cnv_conf: CNV configuration dictionary.
        select_algo: Selection algorithm dictionary.
        driver_genes: Driver genes dictionary.
        can_ids: Canonical IDs dictionary (currently not used for CNV filtering).
        initial_driver_counter: Initial driver counter from SNV processing.
    Returns:
        Tuple of updated driver counter and CNV counter.
    """
    driver_counter = initial_driver_counter
    cnv_counter = 0
    if cnv_file is not None:
        out_path = outputcnv if outputcnv else re.sub(r'\.txt$|\.tsv', DEFAULT_CNV_SUFFIX, os.path.basename(cnv_file))
        logger.info(f"Processing CNV file '{cnv_file}' ...")

        with open(cnv_file, "r", encoding="utf8") as tsvfile, open(out_path, mode='w', newline='') as outtsv:
            tsv_reader = csv.reader(tsvfile, delimiter="\t")
            tsv_writer = csv.writer(outtsv)
            header = ["chrom","loc.start","loc.end","ID","CNt","Geno","logratio","ploidy","call","LOH","gene","driver_status","gene_type"]
            tsv_writer.writerow(header)

            next(tsv_reader, None)

            for row in tsv_reader:
                cnv_counter += 1
                driver_status = 'NotDriver'
                gene_id = row[cnv_conf['gene_id']]
                var_class = row[cnv_conf['class']]
                start = int(row[cnv_conf['start']])
                end = int(row[cnv_conf['end']])
                seglen = end - start

                gene_type = get_gene_type(gene_id, driver_genes)

                # Rule 1: gene-specific rule takes priority
                if (gene_id, 'cnv') in select_algo:
                    is_driver_variant = is_driver_cnv(var_class, seglen, select_algo[(gene_id, 'cnv')])
                    if is_driver_variant:
                        driver_status = 'Driver'

                # Rule 2: gene-type rule (oncogene / tsg / both / unknown)
                if gene_type is not None:
                    key = (gene_type, 'cnv')
                    if key in select_algo:
                        is_driver_variant = is_driver_cnv(var_class, seglen, select_algo[key])
                        if is_driver_variant:
                            driver_status = 'Driver'

                ######################
                # Export
                ######################

                if (driver_status == 'Driver') or (driver_status == 'NotDriver' and debug):
                    if driver_status == 'Driver':
                        driver_counter += 1
                    row.append(driver_status)
                    row.append(gene_type)
                    tsv_writer.writerow(row)
        logger.info(f"CNV processing complete: {cnv_counter} variants processed, {driver_counter - initial_driver_counter} CNV drivers identified")
    return driver_counter, cnv_counter

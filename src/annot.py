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
import yaml
import csv
import gzip

def load_cancer_genes_list(infile: str, format: str = "oncoKB") -> dict:
    """
    Loads the cancer gene list from a TSV file.
    Args:
        infile: Path to the TSV file.
        format: Format of the file (default: "oncoKB").
    Returns:
        Dictionary with keys 'oncogene', 'tsg', 'both', 'unknown' containing lists of gene names.
    """
    oncogene=[]
    tsg=[]
    both=[]
    unknown=[]
    if format == "oncoKB":
        with open(infile, "r", encoding="utf8") as stream:
            tsv_reader = csv.DictReader(stream, delimiter="\t")
            for genes in tsv_reader:
                name = genes["Hugo Symbol"]
                if genes["Is Oncogene"]=="Yes" and genes["Is Tumor Suppressor Gene"]=="No":
                    oncogene.append(name)
                elif genes["Is Oncogene"]=="No" and genes["Is Tumor Suppressor Gene"]=="Yes":
                    tsg.append(name)
                elif genes["Is Oncogene"]=="Yes" and genes["Is Tumor Suppressor Gene"]=="Yes":
                    both.append(name)
                elif genes["Is Oncogene"]=="No" and genes["Is Tumor Suppressor Gene"]=="No":
                    unknown.append(name)
    return {'oncogene': oncogene, 'tsg': tsg, 'both': both, 'unknown': unknown}


def clean_version_number(s: str) -> str:
    """
    Removes version number from a string (e.g., NM or ENSEMBL IDs).
    Args:
        s: Input string.
    Returns:
        String with version number removed.
    """
    return re.sub('\\..*$', '', s)


def get_gene_type(gene_id: str, cancer_atlas: dict) -> str:
    """
    Returns the gene type ('oncogene', 'tsg', 'both', or None) based on the cancer atlas.
    Args:
        gene_id: Gene identifier.
        cancer_atlas: Dictionary from load_cancer_genes_list.
    Returns:
        Gene type or None if not found.
    """
    if gene_id in cancer_atlas['oncogene']:
        return "oncogene"
    elif gene_id in cancer_atlas['tsg']:
        return "tsg"
    elif gene_id in cancer_atlas['both']:
        return "both"
    elif gene_id in cancer_atlas['unknown']:
        return "unknown"
    else:
        return None


def load_canonicals(infile: str, format: str="gtf", clean_ids: bool=False) -> dict:
    """
    Load canonical transcripts from GTF file
    Tested on MANE GTF file
    Of note, some genes can have multiple canonical transcripts (Mane SELECT and CLINICAL_PLUS)
    Args:
        infile: Path to the GTF file.
        format: Format of the file (default: "gtf").
        clean_ids: Whether to clean version numbers from transcript IDs.
    Returns:
        Dictionary mapping gene IDs to lists of canonical transcript IDs.
    """
    ids_table={}
    black_list=[]
    if format == "gtf":
        if infile.endswith('.gz'):
            stream = gzip.open(infile, 'rt')
        else:
            stream = open(infile, mode='r')

        for line in stream:
            if not line.startswith("#"):
                cols = line.strip().split("\t")
                s = next(csv.reader([cols[8]], delimiter=' '))
                tids = s[s.index('transcript_id') + 1].rstrip(";")
                if clean_ids:
                    tids=clean_version_number(tids)
                gids = s[s.index('gene_id') + 1].rstrip(";")

                ## Could have several Ids from MANE Select and Clinical Plus annotation
                if gids in ids_table:
                    if tids not in ids_table[gids]:
                        ids_table[gids].append(tids)
                elif tids != "":
                    ids_table[gids] = [tids]

    return ids_table




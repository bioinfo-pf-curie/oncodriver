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
DEFAULT_VCF_SUFFIX = '_oncoDriver.vcf.gz'

import cyvcf2
import gzip
import logging
import re
import numpy as np
import os
import sys

from .annot import clean_version_number, get_gene_type

def is_driver(variant: dict, db_info: dict, conf_select: list, conf_vcf: dict) -> list:
    """
    Determine if a variant is a driver mutation.
    Args:
        variant: Variant annotation dictionary.
        db_info: Variant INFO fields.
        conf_select: Selection rules from config.
        conf_vcf: VCF configuration.
    Returns:
        List [bool, str or None]: [is_driver, decision_class]
    """
    ## Separate information using sep
    var_info = variant[conf_vcf['annot_info']]
    if conf_vcf['sep'] is not None:
        var_info = var_info.split(conf_vcf['sep'])

    # Hotspot
    in_databases = is_hotspot(variant, infos=db_info, flags=conf_vcf['cancer_db'])
    for r in conf_select:
        if 'is_hotspot' in r and r['is_hotspot'] is True and not in_databases:
            continue
        for v in var_info:
            if v in r['var_classes']:
                return [True, v]
    return [False, None]


def is_polym(v, infos, flags:dict, val: float) -> bool:
    """
    Check if a variant is in a genomeDb with a MAF > val
    Args:
        v: Variant object.
        infos: INFO field dictionary.
        flags: Dictionary of genome databases to check.
        val: MAF threshold. 
    Returns:
        True if variant is polymorphic, False otherwise.
    """    

    fdb = []
    for db in flags:
        for x in flags[db]:
            fdb.append(x)
    
    sub_INFO = subset_INFO(infos, keys=fdb)
    for key in sub_INFO:
        if type(sub_INFO[key]) is tuple:
            for i in sub_INFO[key]:
                if i is not None and i != ".":
                    if float(i) >= float(val):
                        return True
        elif sub_INFO[key] is not None and sub_INFO[key] != "." and sub_INFO[key] != "NA":
            if float(sub_INFO[key]) >= float(val):
                return True
    return False


def is_hotspot(v, infos: dict, flags: list) -> bool:
    """
    Check if a variant is annotated as a cancer hotspot.
    Args:
        v: Variant object.
        infos: INFO field dictionary.
        flags: List of hotspot flags to check.
    Returns:
        True if hotspot, False otherwise.
    """
    subINFO = subset_INFO(infos, keys=flags)
    for key in subINFO:
        if subINFO[key] is not None and subINFO[key] != ".":
            return True
    return(False)


def subset_INFO(annot, keys: list):
    """
    Subset the annotation information to a few key values
    Args:
        annot: Annotation dictionary or list of dictionaries.
        keys: List of keys to subset.
    Returns:
        Subset of annotation information.
    """
    
    if isinstance(annot, list):
        subset_info = []
        for i in range(0, len(annot)):
            z = dict((k, annot[i][k]) for k in keys if k in annot[i])
            if len(z) > 0:
                subset_info.append(z)
    else:
        subset_info = dict((k, annot[k]) for k in keys if k in annot)
    return(subset_info)


def get_INFO(INFO: str) -> list:
    """
    Format the INFO field from snpEff and return a list of dict.
    Args:
        INFO: INFO field string.
    Returns:
        List of dictionaries with annotation info.
    """
    if INFO is not None:
        annotTag = INFO.split(',')
        annotInfo = []
        for i in range(0, len(annotTag)):
            annot = annotTag[i].split('|')
            dictannot = {i: annot[i] for i in range(0, len(annot))}
            annotInfo.append(dictannot)
        return(annotInfo)


def get_tag(v, tag: str) -> np.ndarray:
    """
    Get a tag value from either the format field or the info field.
    Return a 2D numpy array.
    Args:
        v: Variant object.
        tag: Tag name.
    Returns:
        Numpy array of tag values.
    """
    # First check in FORMAT field
    if tag in v.FORMAT:
        val = v.format(tag)

    # Otherwise, check in INFO field
    if tag not in v.FORMAT or val is None:
        val = v.INFO.get(tag)

    if type(val) != np.ndarray:
        val = np.array([val], float)
    else:
        val = val.astype('float')

    return(val)


def process_vcf(vcf_file: str, sample: str, output: str, debug: bool, min_vaf: float, min_depth: int, min_alt_depth: int, max_maf: float, use_canonical: bool, strict: bool,
                vcf_conf: dict, select_algo: dict, driver_genes: dict, can_ids: dict) -> tuple[int, int]:
    """
    Process VCF file and return driver and SNV counters.
    Args:
        vcf_file: Path to the VCF file.
        sample: Sample name to focus on.
        output: Output file name.
        debug: Debug mode flag.
        min_vaf: Minimum variant allele frequency.
        min_depth: Minimum sequencing depth.
        min_alt_depth: Minimum alternative allele depth.
        max_maf: Maximum minor allele frequency.
        use_canonical: Use canonical transcripts flag.
        strict: Strict mode flag.
        vcf_conf: VCF configuration dictionary.
        select_algo: Selection algorithm dictionary.
        driver_genes: Driver genes dictionary.
        can_ids: Canonical IDs dictionary.
    Returns:
        Tuple of driver counter and SNV counter.
    """
    driver_counter = 0
    snv_counter = 0
    if vcf_file is not None:
        if sample is not None:
            vcf = cyvcf2.VCF(vcf_file)
            count = 0
            for s in vcf.samples:
                count = count + 1
                if str(s) == str(sample):
                    vcf = cyvcf2.VCF(vcf_file, samples=sample)
                elif count == len(vcf.samples):
                    print("Error: Name of the sample incorrect\n")
                    sys.exit(-1)
        else:
            vcf = cyvcf2.VCF(vcf_file)

        # Sample name
        if len(vcf.samples) > 1:
            sys.stderr.write("Error: " + str(len(vcf.samples)) +
                             " sample detected. This version is designed for a single sample ! Use --sample argument.\n")
            sys.exit(-1)
        
        # Outputs
        vcf.add_info_to_header({'ID': 'ONCODRIVER', 'Description': 'Oncogenic driver decision',
                                'Type': 'Character', 'Number': '1'})

        if output:
            wx = cyvcf2.Writer(output, vcf)
        else:
            wx = cyvcf2.Writer(re.sub(r'\.vcf$|\.vcf.gz$|\.bcf', DEFAULT_VCF_SUFFIX, os.path.basename(vcf_file)), vcf)
        
        for variant in vcf:
            __DEBUGINFO__ = None
            __DRIVERINFO__ = None
            DECISION = None
            snv_counter += 1

            if snv_counter % 1000 == 0:
                logging.info(f"## {snv_counter}")
 
            try:
                db_info = dict(variant.INFO)

                # Get annotation INFO as a list of dict
                annot_info = get_INFO(variant.INFO.get(vcf_conf['tag']))
                if annot_info is None:
                    raise Exception("Annotation tag [" + vcf_conf['tag'] + "] not found !")
                
                ##############################
                # FORMAT - technical filters
                ##############################
                tech_filtered=False
                
                # Variant Allele Frequency
                if min_vaf is not None:
                    fval = get_tag(variant, vcf_conf['vaf'])
                    if fval is not None and len(fval[fval <= min_vaf]) == len(variant.ALT):
                        __DEBUGINFO__ = ",".join(x for x in [__DEBUGINFO__, "VAF"] if x is not None)
                        tech_filtered = True
                        if not debug:
                            continue
                    elif fval is None:
                        raise Exception("VAF tag [" + vcf_conf['vaf'] + "] not found !")

                # Sequencing Depth
                if min_depth is not None:
                    dval = get_tag(variant, vcf_conf['depth'])
                    if dval is not None and len(dval[dval <= min_depth]) == len(variant.ALT):
                        __DEBUGINFO__ = ",".join(x for x in [__DEBUGINFO__, "DEPTH"] if x is not None)
                        tech_filtered = True
                        if not debug:
                            continue
                    elif dval is None:
                        raise Exception("DEPTH tag [" + vcf_conf['depth'] + "] not found !")

                # Alternative allele Depth
                if min_alt_depth is not None:
                    ad = get_tag(variant, vcf_conf['alt_depth'])
                    # case where AD = REF + ALTs
                    if len(ad) == (len(variant.ALT) + 1):
                        ad = ad[1:]
                        
                    if ad is not None and len(ad[ad <= min_alt_depth]) == len(variant.ALT):
                        __DEBUGINFO__ = ",".join(x for x in [__DEBUGINFO__, "ALTDEPTH"] if x is not None)
                        tech_filtered = True
                        if not debug:
                            continue
                    elif ad is None:
                        raise Exception("ALT DEPTH tag [" + vcf_conf['alt_depth'] + "] not found !")

                 # Polymorphisms
                if max_maf is not None:
                    if is_polym(variant, infos=db_info, flags=vcf_conf['polym_db'], val=max_maf):
                        __DEBUGINFO__ = ",".join(x for x in [__DEBUGINFO__, "POLYM"] if x is not None)
                        tech_filtered = True
                        if not debug:
                            continue
                    
                ######################
                # Annotation INFO
                ######################

                if not tech_filtered :
                    for annot in annot_info:
                        try:
                            ## Extract info from variant
                            gene_id = annot[vcf_conf['gene_id']]
                            transcript_id=clean_version_number(annot[vcf_conf['transcript_id']])
                        except KeyError as e:
                            logging.warning(f"Annotation key {e} not found, skipping annotation")
                            continue
                        gene_type = get_gene_type(gene_id, driver_genes)
                        is_driver_variant = False

                        ## Genes specified in the decision algorithm
                        if (gene_id, 'snv') in select_algo:
                            if (use_canonical and (can_ids.get(gene_id) is not None) and transcript_id not in can_ids[gene_id]):
                                if transcript_id not in can_ids[gene_id]:
                                    __DEBUGINFO__ = ",".join(x for x in [__DEBUGINFO__, "NON_CANONICAL"] if x is not None)
                                continue
                            elif not use_canonical or (use_canonical and (can_ids.get(gene_id) is not None) and transcript_id in can_ids[gene_id]):
                                 is_driver_variant, decision = is_driver(annot, db_info, select_algo[(gene_id, 'snv')], vcf_conf)
                                 __DRIVERINFO__ = ",".join(x for x in [__DRIVERINFO__, transcript_id, gene_type, decision ] if x is not None)
                                 if not is_driver_variant:
                                    if debug:
                                        __DEBUGINFO__ = ",".join(x for x in [__DEBUGINFO__, "NON_DRIVER"] if x is not None)
                                        continue
                                 else:
                                    if debug:
                                        __DEBUGINFO__ = str(is_driver_variant)
                                    driver_counter += 1
                                    break

                        ## Driver genes from cancer gene list
                        elif gene_type is not None:
                            if (use_canonical and (can_ids.get(gene_id) is not None) and transcript_id not in can_ids[gene_id]):
                                if transcript_id not in can_ids[gene_id]: 
                                    __DEBUGINFO__ = ",".join(x for x in [__DEBUGINFO__, "NON_CANONICAL"] if x is not None)
                                continue
                            elif not use_canonical or (use_canonical and (can_ids.get(gene_id) is not None) and transcript_id in can_ids[gene_id]):
                                is_driver_variant, decision = is_driver(annot, db_info, select_algo[(gene_type, 'snv')], vcf_conf)
                                __DRIVERINFO__ = ",".join(x for x in [__DRIVERINFO__, transcript_id, gene_type, decision ] if x is not None)
                                if not is_driver_variant:
                                    if debug:
                                        __DEBUGINFO__ = ",".join(x for x in [__DEBUGINFO__, "NON_DRIVER"] if x is not None)
                                        continue
                                else:
                                    if debug:
                                        __DEBUGINFO__ = str(is_driver_variant)
                                    driver_counter += 1
                                    break
                        elif debug:
                            __DEBUGINFO__ = ",".join(x for x in [__DEBUGINFO__, "NON_CANCER_GENES"] if x is not None)
                            break

                if debug:
                    if __DEBUGINFO__ is not None:
                        variant.INFO['ONCODRIVER'] = __DEBUGINFO__
                    wx.write_record(variant)
                elif is_driver_variant is True:
                    variant.INFO['ONCODRIVER'] = __DRIVERINFO__
                    wx.write_record(variant)

            except Exception as e:
                warnflag = str(variant.CHROM) + ":" + str(variant.start+1) + "[" + str(variant.REF) + "/" + str(variant.ALT[0]) + "]"
                logging.warning(f"Warning : variant {warnflag} raises an error: \"{e}\". Skipped ...")
                if strict:
                    logging.error("Error - strict mode activated - exit")
                    sys.exit(-1)

        wx.close()
        vcf.close()
    return driver_counter, snv_counter



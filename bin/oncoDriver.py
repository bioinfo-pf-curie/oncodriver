#!/usr/bin/env python
# coding: utf-8
#
#  This file is part of oncoDriver software.
#
#  Copyright (c) 2024 - Institut Curie
#
#  File author(s):
#      Nicolas Servant <nicolas.servant@curie.fr>
#
#  Distributed under the terms of the CeCILL-B license.
#  The full license is in the LICENSE file, distributed with this software.
#
##############################################################################

__version__ = '1.0.0dev'

"""
This script is designed to apply a decision tree on a VCF file
and to select driver mutations from a VCF file.

python oncoDriver.py -i ${VCF} \
  --config ${DECISION_TREE} \
  --use_canonical \
  --canonical_reference ${TABLE}  > oncoDriver_results.log
"""

import argparse
import sys
import warnings
import re
import yaml
import numpy as np
import os.path
import csv
from datetime import date
import cyvcf2
import gzip

"""
Load the decision rules and return a hash table with (gene_type,var_type)
"""

def load_config(infile):
    cfg={}
    with open(infile, 'r') as stream:
        try:
            data=yaml.safe_load(stream)
            cfg['vcf'] = data['vcf'] 
            d={}
            for info in data['select']:
                k=(info['gene_type'],info['var_type'])
                if ( k in d ):
                    d[k].append(info)
                else:
                    d[k]=[info]
            cfg['select'] = d
            return cfg
        except:
            raise

"""
Load cancer gene list
Supported format : oncoKB
"""

def load_cancer_genes_list(infile, format="oncoKB"):

    oncogene=[]
    tsg=[]
    both=[]
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
    return {'oncogene': oncogene, 'tsg': tsg, 'both': both}


"""
Remove version number from a string (NM or ENSEMBL Ids)
"""

def clean_version_number(s):
    return re.sub('\.[0-9]*$', '', s)

"""
Combined cancer gene atlas with gene id
and retur the gene type 'oncogene'/'tsg'/'both'
"""

def get_gene_type(gene_id, cancer_atlas):

    if gene_id in cancer_atlas['oncogene']:
        return "oncogene"
    elif gene_id in cancer_atlas['tsg']:
        return "tsg"
    elif gene_id in cancer_atlas['both']:
        return "both"
    else:
        return None

"""
Load canonical transcripts from GTF file
Tested on MANE GTF file
Of note, some genes can have multiple canonical transcripts (Mane SELECT and CLINICAL_PLUS)
"""

def load_canonicals(infile, format="gtf", clean_ids=False):

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

"""
Apply the decision tree to a variant
return is the variant is driver or not
"""

def is_driver(variant, db_info, conf_select, conf_vcf):

    ## TODO use SEP ? splice_region_variant&non_coding_transcript_exon_variant

    # Hotspot
    in_databases = is_hotspot(variant, infos=db_info, flags=conf_vcf['cancer_db'])
    for r in conf_select:
        if 'is_hotspot' in r and not in_databases:
            return False
        if variant[conf_vcf['annot_info']] in r['var_classes']:
            return True
    return False


"""
Check if a variant is in a genomeDb with a MAF > val
"""

def is_polym(v, infos, flags, val):

    ## TODO - add check on gnomAD=True ?
    
    ## Flatten polym info
    fdb = []
    for db in flags:
        for x in flags[db]:
            fdb.append(x)
    
    sub_INFO = subset_INFO(infos, keys=fdb)
    for key in sub_INFO:
        if type(sub_INFO[key]) is tuple:
            for i in sub_INFO[key]:
                print(sub_INFO[key])
                if i is not None and i != ".":
                    if float(i) >= float(val):
                        return True
        elif sub_INFO[key] is not None and sub_INFO[key] != "." and sub_INFO[key] != "NA":
            if float(sub_INFO[key]) >= float(val):
                return True
    return False


"""
Check if a variant is annotated as a cancer hotspot
"""

def is_hotspot(v, infos, flags):

    subINFO = subset_INFO(infos, keys=flags)
    for key in subINFO:
        if subINFO[key] is not None and subINFO[key] != ".":
            return True
    return(False)


"""
Subset the annotation information to a few key values
"""

def subset_INFO(annot, keys):
    
    if isinstance(annot, list):
        subset_info = []
        for i in range(0, len(annot)):
            z = dict((k, annot[i][k]) for k in keys if k in annot[i])
            if len(z) > 0:
                subset_info.append(z)
    else:
        subset_info = dict((k, annot[k]) for k in keys if k in annot)
    return(subset_info)


"""
Format the INFO field from snpEff and return a list of dict
ie. snpEff
"""

def get_INFO(INFO):

    if INFO is not None:
        annotTag = INFO.split(',')
        annotInfo = []
        for i in range(0, len(annotTag)):
            annot = annotTag[i].split('|')
            dictannot = {i: annot[i] for i in range(0, len(annot))}
            annotInfo.append(dictannot)
        return(annotInfo)


"""
Get a tag value from either the format field or the info field
Return a 2D numpy array
"""

def get_tag(v, tag):
    
    # First check in FORMAT field
    if tag in variant.FORMAT:
        val = variant.format(tag)

    # Otherwise, check in INFO field
    if tag not in variant.FORMAT or val is None:
        val = variant.INFO.get(tag)

    if type(val) != np.ndarray:
        val = np.array([val], float)
    else:
        val = val.astype('float')

    return(val)


"""
Parse inputs
"""

def args_parse():

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", "--vcf", required=True, help="Input file (.vcf, .vcf.gz, .bcf)")

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
    parser.add_argument("--oncoKB_genes", help="Cancer gene list from oncoKB", type=str, required=True)
    
    # Others
    parser.add_argument("--verbose", help="Active verbose mode", action="store_true")
    parser.add_argument("--debug", help="Export original VCF with TMB_FILTER tag", action="store_true")
    parser.add_argument("--version", help="Version number", action='version', version="%(prog)s ("+__version__+")")

    args = parser.parse_args()
    return (args)


if __name__ == "__main__":

    args = args_parse()

    # Load Data
    if args.sample is not None:
        vcf = cyvcf2.VCF(args.vcf)
        count = 0
        for sample in vcf.samples:
            count = count + 1
            if str(sample) == str(args.sample):
                vcf = cyvcf2.VCF(args.vcf, samples=args.sample)
            elif count == len(vcf.samples):
                print("Error: Name of the sample incorrect\n")
                sys.exit(-1)
    else:
        vcf = cyvcf2.VCF(args.vcf)

    # Sample name
    if len(vcf.samples) > 1:
        sys.stderr.write("Error: " + str(len(vcf.samples)) +
                         " sample detected. This version is designed for a single sample ! Use --sample argument.\n")
        sys.exit(-1)
        
    # Outputs
    if args.debug:
        vcf.add_info_to_header({'ID': 'ONCODRIVER', 'Description': 'Is an oncogenic driver',
                                'Type': 'Character', 'Number': '1'})

    wx = cyvcf2.Writer(re.sub(r'\.vcf$|\.vcf.gz$|\.bcf',
                              '_oncoDriver.vcf.gz', os.path.basename(args.vcf)), vcf)
    
    # Load config
    conf = load_config(args.config)
    select_algo = conf['select']
    vcf_conf = conf['vcf']
 
    # Load cancer genes
    driver_genes = load_cancer_genes_list(args.oncoKB_genes)

    # Load canonical Ids
    can_ids = None
    if args.use_canonical:
        can_ids = load_canonicals(args.canonical_ids, clean_ids=True)
        
    # Init counter
    driver_counter = 0
    var_counter = 0
        
    for variant in vcf:
        __DEBUGINFO__ = None
        var_counter += 1
        if (var_counter % 1000 == 0 and args.verbose):
            print ("## ", var_counter)
 
        try:
            db_info = dict(variant.INFO)

            # Get annotation INFO as a list of dict
            annot_info = get_INFO(variant.INFO.get(vcf_conf['tag']))
            
            ##############################
            # FORMAT - technical filters
            ##############################
            tech_filtered=False
            
            # Variant Allele Frequency
            if args.min_vaf is not None:
                fval = get_tag(variant, vcf_conf['vaf'])
                if fval is not None and len(fval[fval <= args.min_vaf]) == len(variant.ALT):
                    __DEBUGINFO__ = ",".join(x for x in [__DEBUGINFO__, "VAF"] if x is not None)
                    tech_filtered = True
                    if not args.debug:
                        continue

            # Sequencing Depth
            if args.min_depth is not None:
                dval = get_tag(variant, vcf_conf['depth'])
                if dval is not None and len(dval[dval <= args.min_depth]) == len(variant.ALT):
                    __DEBUGINFO__ = ",".join(x for x in [__DEBUGINFO__, "DEPTH"] if x is not None)
                    tech_filtered = True
                    if not args.debug:
                        continue

            # Alternative allele Depth
            if args.min_alt_depth is not None:
                ad = get_tag(variant, vcf_conf['alt_depth'])
                # case where AD = REF + ALTs
                if len(ad) == (len(variant.ALT) + 1):
                    ad = ad[1:]
                    
                if ad is not None and len(ad[ad <= args.min_alt_depth]) == len(variant.ALT):
                    __DEBUGINFO__ = ",".join(x for x in [__DEBUGINFO__, "ALTDEPTH"] if x is not None)
                    tech_filtered = True
                    if not args.debug:
                        continue

             # Polymorphisms
            if args.max_maf is not None:
                if is_polym(variant, infos=db_info, flags=vcf_conf['polym_db'], val=args.max_maf):
                    __DEBUGINFO__ = ",".join(x for x in [__DEBUGINFO__, "POLYM"] if x is not None)
                    tech_filtered = True
                    if not args.debug:
                        continue
                    
            ######################
            # Annotation INFO
            ######################

            if not tech_filtered :
                for annot in annot_info:
                    gene_id = annot[vcf_conf['gene_id']]
                    transcript_id=clean_version_number(annot[vcf_conf['transcript_id']])
                
                    gene_type = get_gene_type(gene_id, driver_genes)
                    is_driver_mut = False
                    if gene_type is not None:
                        if (args.use_canonical and transcript_id not in can_ids[gene_id]):
                            __DEBUGINFO__ = ",".join(x for x in [__DEBUGINFO__, "NON_CANONICAL"] if x is not None)
                            continue
                        elif not args.use_canonical or (args.use_canonical and transcript_id in can_ids[gene_id]):
                            is_driver_mut = is_driver(annot, db_info, select_algo[(gene_type, 'snv')], vcf_conf)
                            if not is_driver_mut:
                                if args.debug:
                                    __DEBUGINFO__ = ",".join(x for x in [__DEBUGINFO__, "NON_DRIVER"] if x is not None)
                                    continue
                            else:
                                driver_counter += 1
                                break
                    elif args.debug:
                        __DEBUGINFO__ = ",".join(x for x in [__DEBUGINFO__, "NON_CANCER_GENES"] if x is not None)
                        break

            ######################
            # Export
            ######################

            if args.debug:
                if __DEBUGINFO__ is not None:
                    variant.INFO["ONCODRIVER"] = __DEBUGINFO__
                else:
                    variant.INFO["ONCODRIVER"] = str(is_driver_mut)
                print(variant)
                wx.write_record(variant)
            elif is_driver_mut is True:
                print(variant)
                wx.write_record(variant)

        except:
            warnflag = str(variant.CHROM) + ":" + str(variant.start) + "-" + str(variant.end)
            warnings.warn("Warning : variant {} raises an error. Skipped so far ...".format(warnflag))
            raise

    wx.close()
    vcf.close()

    # Output
    print("## oncoDriver version=", __version__, file=sys.stderr)
    print("## When=", date.today(), file=sys.stderr)
    print("")
    print("## Input=", args.vcf, file=sys.stderr)
    print("## Sample=", args.sample, file=sys.stderr)
    print("")
    print("## Config=", args.config, file=sys.stderr)
    print("## oncoKB genes=", args.oncoKB_genes, file=sys.stderr)
    print("## Use canonical=", args.use_canonical, file=sys.stderr)
    print("## Canonical genes=", args.canonical_ids, file=sys.stderr)
    print("")
    print("## Total variants=", var_counter, file=sys.stderr)
    print("## Driver variants=", driver_counter, file=sys.stderr)

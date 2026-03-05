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
import yaml

logger = logging.getLogger(__name__)


def load_configuration(config_file: str) -> dict:
    """
    Loads the decision rules from a YAML configuration file.
    Args:
        config_file: Path to the YAML configuration file.
    Returns:
        Configuration dictionary with 'vcf', 'cnv', and 'select' keys.
    Raises:
        ValueError: If the YAML file cannot be read or parsed.
    """
    cfg = {}
    with open(config_file, 'r') as stream:
        try:
            data = yaml.safe_load(stream)
            for key in ('vcf', 'cnv', 'select'):
                if key not in data:
                    raise ValueError(f"Missing required section '{key}' in config file: {config_file}")
            cfg['vcf'] = data['vcf']
            cfg['cnv'] = data['cnv']
            d={}
            for info in data['select']:
                if 'gene_id' in info:
                    for gid in info['gene_id'].split("|"):
                        k=(gid.strip(), info['var_type'])
                        gid_info = info.copy()
                        gid_info['gene_id'] = gid.strip()
                        if (k in d):
                            d[k].append(gid_info)
                        else:
                            d[k]=[gid_info]
                elif 'gene_type' in info:
                    for gt in info['gene_type'].split("|"):
                        k=(gt,info['var_type'])
                        gt_info = info.copy()
                        gt_info['gene_type']=gt
                        if ( k in d ):
                            d[k].append(gt_info)
                        else:
                            d[k]=[gt_info]
            cfg['select'] = d
            return cfg
        except yaml.YAMLError as e:
            raise ValueError(f"Error - unable to read {config_file}: {e}")


def print_config(cfg: dict) -> None:
    """
    Prints the decision algorithm configuration in a human-readable manner.
    Args:
        cfg: Configuration dictionary to print.
    """
    for section, val in cfg.items():
        if isinstance(val, dict):
            logger.info(f"## {section}->")
            for j, val2 in val.items():
                if isinstance(val2, list):
                    for item in val2:
                        logger.info(f"## - {j} = {item}")
                else:
                    logger.info(f"## - {j} = {val2}")
        else:
            logger.info(f"## {section} = {val}")

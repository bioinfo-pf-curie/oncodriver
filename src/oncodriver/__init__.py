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

__version__ = "1.3.0"

from .annot import (
    load_cancer_genes_list,
    load_canonicals,
    clean_version_number,
    get_gene_type,
)
from .utils import load_configuration, print_config
from .snv_process import process_vcf
from .cnv_process import process_cnv

__all__ = [
    "__version__",
    "load_cancer_genes_list",
    "load_canonicals",
    "clean_version_number",
    "get_gene_type",
    "load_configuration",
    "print_config",
    "process_vcf",
    "process_cnv",
]

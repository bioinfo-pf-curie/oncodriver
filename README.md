# oncoDriver

[![Install with conda](https://img.shields.io/badge/install%20with-conda-brightgreen.svg)](https://conda.anaconda.org/anaconda)

This package is designed to filter and select driver alterations for cancer genes from both SNVs (vcf) and CNVs (tsv).
Variants can be filtered following both technical and biological rules defined in a configuration file.

## Installation

### From source (pip)

Python ≥ 3.10 is required.

```bash
git clone https://gitlab.com/nservant/oncodriver.git
cd oncodriver
pip install .
```

To also install the optional testing and linting dependencies:

```bash
pip install ".[test,lint]"
```

### With conda

A ready-to-use conda environment is provided via [`environment.yml`](environment.yml):

```bash
conda env create -f environment.yml
conda activate oncodriver
pip install .
```

The conda environment pins the key bioinformatics dependencies (`cyvcf2`, `pyyaml`, `mosdepth`, `pybedtools`) to tested versions.

---

## Testing

The test suite uses [pytest](https://docs.pytest.org/) with [pytest-cov](https://pytest-cov.readthedocs.io/) for coverage reporting.

### Install test dependencies

```bash
pip install ".[test]"
```

### Run all tests

```bash
pytest
```

By default (see [`pyproject.toml`](pyproject.toml)), pytest:

- discovers tests in the `tests/` directory
- measures coverage of the `oncodriver` package (`--cov=oncodriver`)
- prints a term-missing coverage report (`--cov-report=term-missing`)

### Run a specific test module

```bash
pytest tests/test_snv_process.py
pytest tests/test_cnv_process.py
pytest tests/test_cli.py
```

### Run with verbose output

```bash
pytest -v
```

### Generate an HTML coverage report

```bash
pytest --cov=oncodriver --cov-report=html
# then open htmlcov/index.html
```

---

## Usage

```
usage: oncodriver [-h] [-i VCF] [--cnv CNV] --config CONFIG [--sample SAMPLE] [--min_vaf MIN_VAF] [--max_maf MAX_MAF] [--min_depth MIN_DEPTH] [--min_alt_depth MIN_ALT_DEPTH] [--use_canonical] [--canonical_ids CANONICAL_IDS] --driver_genes DRIVER_GENES [--output OUTPUT] [--outputcnv OUTPUTCNV] [--verbose] [--debug] [--strict] [--version]

options:
  -h, --help                    show this help message and exit
  -i VCF, --vcf VCF             Input file (.vcf, .vcf.gz, .bcf) (default: None)
  --cnv CNV                     Input file for CNV (segments.transformed_annot_oncokb.txt) (default: None)
  --config CONFIG               Decision tree config file (default: None)
  --sample SAMPLE               Specify the sample ID to focus on (default: None)
  --min_vaf MIN_VAF             Filter variants with Allelic Ratio <= vaf (default: None)
  --max_maf MAX_MAF             Filter variants with MAF > maf (default: None)
  --min_depth MIN_DEPTH         Filter variants with depth < min_depth (default: None)
  --min_alt_depth MIN_ALT_DEPTH Filter variants with alternative allele depth <= min_alt_depth (default: None)
  --use_canonical               Use canonical transcript only (default: False)
  --canonical_ids CANONICAL_IDS File with canonical IDs information (default: None)
  --driver_genes DRIVER_GENES   Cancer gene list from oncoKB (default: None)
  --output OUTPUT               Output file name (default: None)
  --outputcnv OUTPUTCNV         CNV output file name (default: None)
  --verbose                     Active verbose mode (default: False)
  --debug                       Export original VCF with ONCODRIVER tag (default: False)
  --strict                      Raise an error instead of warnings (default: False)
  --version                     Version number
																																							
```

### Command line

```
oncodriver \
  --vcf D1326R01_vs_D1326R02_Mutect2_filtered_pass_norm_COSMIC_ICGC_CancerHotspots_GnomAD_dbNSFP.vcf.gz \
  --sample D1326R02 \
  --config ./config/pathogenic_variants.yml \
  --oncoKB_genes ./assets/cancerGeneList.tsv \
  --canonical_ids MANE.GRCh38.v1.3.refseq_genomic.gtf.gz \
  --use_canonical
```

### Inputs

SNVs variants are expected to be in `vcf.gz` format. The variants must be annotated with `snpEff` or `VEP`.

The CNV variants are expected to be in a tsv format with the following information ;  
"chrom","loc.start","loc.end","ID","CNt","Geno","logratio","ploidy","call","LOH","gene"

This file can be easily generated from `Facets` outputs.

### Configuration file

The decision algorithm is written in a `yaml` file which can be customized according to your criteria.
A few section are mandatory.

- `vcf`: define where to find some information in the `vcf` file format

Here is an example to parse the annotation from `snpEff`:

```
- vcf
  tag: 'ANN'
  sep: '&'
  vaf: 'AF'
 depth: 'DP'
  alt_depth: 'AD'
  annot_info: 1
  gene_id: 3
  transcript_id : 6
  cancer_db:
    - 'ICGC'
    - 'CancerHotspots'
    - 'COSMIC'
  polym_db:
    gnomAD:
      - AF
```

- `cnv`: define how to parse the txt file with copy number information

```
cnv:
  gene_id: 10
  class: 8
  start: 1
  end: 2
```

- `select`: describe the rules to consider a variant as driver.

Those rules usually contain the following information:
- `var_type`: must be `snv` or `cnv`
- `gene_type`: must be `oncogenes`, `tsg`, `both` or `unknown`. Of note, multiple gene type can be considered, separated by a `|` character.
- `var_classes`: list here all type of variants to consider as driver
- `is_hotspot`: define whether the variant must be (or not) annotated at least in of the databases define in `cancer_db`


For instance ; *I would like to consider as driver all misense variants occuring on a oncogene, and known in the cancer databases*
will be translated into ;

```
  - var_type: 'snv'
    gene_type: 'oncogene'
    var_classes:
      - missense_variant
    is_hotspot: true
```

In order to rescue a few variants which are on specific genes, such as for instance TERT mutation in the promoter region, the selection can also be applied by `gene_id`.
For instance, to select the variants occuring in the TERT promoter region, you can use;

```
  - var_type: 'snv'
    gene_id: 'TERT'
    var_classes:
      - upstream_gene_variant
```

### Output

The default output is the list of drivers mutations in vcf/txt format.

### `--debug` mode

In `--debug` mode, all variants are outputed with a flag explaining why it has been selected or filtered.
The following flag are available :

- NON_CANONICAL : The variant is not on the canonical transcript (see `--use_canonical` and `--canonical_ids`)
- NON_CANCER_GENES : The variant is on a gene which is not a driver gene (see `--driver_genes`)
- NON_DRIVER : The variant is not driver as defined in the configuration file
- VAF : The variant allele frequency is lower than expected (see `--min_vaf`)
- DEPTH : The variant depth is not enough covered (see `--min_depth`)
- ALTDEPTH : The alternative allele is not enough covered (see `--min_alt_depth`)
- POLYM : The variant is a polymorphism


## AI Disclosure: Augmented
This project is **AI-augmented** and utilized AI (e.g., Claude) to:
* **Generate** boilerplate code and specific utility functions.
* **Refactor** existing code for better performance and readability.
* **Draft** unit tests and technical documentation.
**Verification:** Every AI-generated contribution was manually reviewed, debugged, and integrated into the final codebase.


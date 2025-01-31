# oncoDriver

[![Install with conda](https://img.shields.io/badge/install%20with-conda-brightgreen.svg)](https://conda.anaconda.org/anaconda)


*Oncodriver*

Filter variants (single nucleotide variants and/or copy number) and select drivers events of clinical interest.

The main idea is *oncodriver* is to apply a selection algorithm describing which variants should be considered as clinically relevant.
The selection rules can be different according to the gene types (oncogenes, tumor suppressor, etc.), and whether the variants is known in cancer databases.
The rules are defined in a `yaml` file, and vary according to the variant annotations.

### Quick Help

```
usage: oncoDriver.py [-h] -i VCF [--cnv CNV] --config CONFIG [--sample SAMPLE] [--min_vaf MIN_VAF] [--max_maf MAX_MAF] [--min_depth MIN_DEPTH] [--min_alt_depth MIN_ALT_DEPTH] [--use_canonical]
                     [--canonical_ids CANONICAL_IDS] --driver_genes DRIVER_GENES [--output OUTPUT] [--verbose] [--debug] [--version] [--outputcnv OUTPUTCNV]

options:
  -h, --help                     show this help message and exit
  -i VCF, --vcf VCF              Input file (.vcf, .vcf.gz, .bcf) (default: None)
  --cnv CNV                      Input file for CNV (segments.transformed_annot_oncokb.txt) (default: None)
  --config CONFIG                Decision tree config file (default: None)
  --sample SAMPLE                Specify the sample ID to focus on (default: None)
  --min_vaf MIN_VAF              Filter variants with Allelic Ratio <= vaf (default: None)
  --max_maf MAX_MAF              Filter variants with MAF > maf (default: None)
  --min_depth MIN_DEPTH          Filter variants with depth < min_depth (default: None)
  --min_alt_depth MIN_ALT_DEPTH  Filter variants with alternative allele depth <= min_alt_depth (default: None)
  --use_canonical                Use canonical transcript only (default: False)
  --canonical_ids CANONICAL_IDS  File with canonical IDs information (default: None)
  --driver_genes DRIVER_GENES    Cancer gene list from oncoKB (default: None)
  --output OUTPUT                Output file name (default: None)
  --verbose                      Active verbose mode (default: False)
  --debug                        Export original VCF with TMB_FILTER tag (default: False)
  --version                      Version number
  --outputcnv OUTPUTCNV          CNV output file name (default: None)
```


### Example usage

```
python oncoDriver.py \
  --vcf D1326R01_vs_D1326R02_Mutect2_filtered_pass_norm_COSMIC_ICGC_CancerHotspots_GnomAD_dbNSFP.vcf.gz \
  --sample D1326R02 \
  --config ./config/pathogenic_variants.yml \
  --oncoKB_genes ./assets/cancerGeneList.tsv \
  --canonical_ids MANE.GRCh38.v1.3.refseq_genomic.gtf.gz \
  --use_canonical
```

### Inputs

SNVs variants are expected to be in `vcf.gz` format. The variants must be annotated with `snpEff` or `VEP`.

The CNV variants are expected to be in a txt format.

### Output

The default output is the list of drivers mutations in vcf/txt format.

### `--debug` mode

In `--debug` mode, all variants are outputed with a flag explaining why it has been selected or filtered.
The following flag are available :

- NON_CANONICAL : The variant is not on the canonical transcript (see `--use_canonical` and `--canonical_ids`)
- NON_CANCER_GENES : The variant is on a genes which is not a driver genes (see `--driver_genes`)
- NON_DRIVER : The variant is not driver as defined in the configuration file
- VAF : The variant allele frequency is lower than expected (see `--min_vaf`)
- DEPTH : The variant depth is not enough covered (see `--min_depth`)
- ALTDEPTH : The alternative allele is not enough covered (see `--min_alt_depth`)
- POLYM : The variant is a polymorphism


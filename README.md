# oncoDriver

Filter and select driver alterations for cancer genes

### Usage

```
usage: oncoDriver.py [-h] [-i VCF] [--cnv CNV] --config CONFIG [--sample SAMPLE] [--min_vaf MIN_VAF] [--max_maf MAX_MAF] [--min_depth MIN_DEPTH] [--min_alt_depth MIN_ALT_DEPTH] [--use_canonical] [--canonical_ids CANONICAL_IDS] --driver_genes DRIVER_GENES [--output OUTPUT] [--outputcnv OUTPUTCNV] [--verbose] [--debug] [--strict] [--version]

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
python oncoDriver.py \
  --vcf D1326R01_vs_D1326R02_Mutect2_filtered_pass_norm_COSMIC_ICGC_CancerHotspots_GnomAD_dbNSFP.vcf.gz \
  --sample D1326R02 \
  --config ./config/pathogenic_variants.yml \
  --oncoKB_genes ./assets/cancerGeneList.tsv \
  --canonical_ids MANE.GRCh38.v1.3.refseq_genomic.gtf.gz \
  --use_canonical
```

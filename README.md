# oncoDriver

Filter and select driver alterations for cancer genes

```
python oncoDriver.py \
  --vcf D1326R01_vs_D1326R02_Mutect2_filtered_pass_norm_COSMIC_ICGC_CancerHotspots_GnomAD_dbNSFP.vcf.gz \
  --sample D1326R02 \
  --config ./config/pathogenic_variants.yml \
  --oncoKB_genes ./assets/cancerGeneList.tsv \
  --canonical_ids MANE.GRCh38.v1.3.refseq_genomic.gtf.gz \
  --use_canonical
```

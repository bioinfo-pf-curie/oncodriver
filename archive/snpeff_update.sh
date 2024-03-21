#!/bin/bash

# COSMIC : LEGACY_ID
# icgc : ID
# CancerHotspot : CH
# dbNSFP : SIFT_pred,Polyphen2_HDIV_pred,CADD_phred
# GnomAD : AF


# cosmic : available from https://cancer.sanger.ac.uk/cosmic/download/ with VCF/CosmicCodingMuts.vcf.gz
# gnomad : https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz & https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz.tbi 

# icgc file : https://dcc.icgc.org/api/v1/download?fn=/release_28/Summary/simple_somatic_mutation.aggregated.vcf.gz
# cancerhotspots : http://download.cbioportal.org/cancerhotspots/cancerhotspots.v2.maf.gz

## Dowload dbNSFP :
# https://snpeff.blob.core.windows.net/databases/dbs/GRCh37/dbNSFP_4.1a/dbNSFP4.1a.txt.gz & https://snpeff.blob.core.windows.net/databases/dbs/GRCh37/dbNSFP_4.1a/dbNSFP4.1a.txt.gz.tbi
# https://snpeff.blob.core.windows.net/databases/dbs/GRCh38/dbNSFP_4.1a/dbNSFP4.1a.txt.gz & https://snpeff.blob.core.windows.net/databases/dbs/GRCh38/dbNSFP_4.1a/dbNSFP4.1a.txt.gz.tbi

## Download recent versions of hg19 & hg38 snpeff db into dedicated directories and unzip them 
# https://snpeff.blob.core.windows.net/databases/v5_0/snpEff_v5_0_hg19.zip
# https://snpeff.blob.core.windows.net/databases/v5_0/snpEff_v5_0_hg38.zip

function usage {
    echo -e "Usage: $(basename $1) [options]"
    echo -e "\n [Options] 
    \t-i/--input_dir : input directory
    \t-a/--assembly_information : either hg19 or hg38 
    \t-c/--conda_bin_dir : conda bin directory
    \t-C/--conda_dir : conda directory output
    \t-s/--snpeff_db : snpeff folder containing recent assembly databases
    \t-o/--output_dir : output directory
    \t-v/--vcf_file : vcf variant file to be analized (can be gzipped)
    \t-n/--vcf_samples_name : name of the vcf file (eg :D777E181_vs_D777E133 )
    \t-h/--help : print help\n"

}

if [ ${#} -eq 0 ];then usage ${0} ; exit 1 ; fi

while getopts "i:a:c:C:s:o:v:n:h" option; do
    case "${option}" in
    	        i) in_dir=$(readlink -m ${OPTARG});;
                a) assembly=${OPTARG};;
                c) conda_bin_dir=$(readlink -m ${OPTARG});;
                C) conda_dir=$(readlink -m ${OPTARG});;
                s) snpeff_db=$(readlink -m ${OPTARG});;
                o) out_dir=$(readlink -m ${OPTARG});;
                v) vcf=$(readlink -m ${OPTARG});;
                n) name=${OPTARG};;
		        h) usage ${0} ; exit
    esac
done

# conda_bin_dir="/data/kdi_prod/.kdi/project_workspace_0/1198/acl/02.00/conda/miniconda3/condabin/"
# conda_dir="/data/kdi_prod/.kdi/project_workspace_0/1198/acl/02.00/conda"

# in_dir="/data/kdi_prod/.kdi/project_workspace_0/1198/acl/02.00/snp_eff/databases"
# out_dir="${in_dir}/snpeff/results"

# assembly="hg38"
# snpeff_db="${in_dir}/snpeff/data"

# vcf="/data/users/egirard1/EXPRESS-21/analysis/WES/Mutect2/D777E181_vs_D777E133_Mutect2_filtered_pass_norm.vcf.gz"
# name="D777E181_vs_D777E133"

# Conda env creation

export PATH=${conda_bin_dir}:$PATH

## Install recent conda (v5.1) of snpeff and snpsift

mkdir -p ${in_dir}/snpeff/ ${in_dir}/snpsift ${conda_dir}/conda_snpeff/  ${conda_dir}/conda_snpsift/ 

if [ ! -d "${conda_dir}/conda_snpeff/snpeff_v5.1" ] ; then mamba env create -f ${in_dir}/snpeff/snpeff.yaml -p ${conda_dir}/conda_snpeff/snpeff_v5.1 ; fi

if [ ! -d "${conda_dir}/conda_snpsift/snpsift_v5.1" ] ; then mamba env create -f ${in_dir}/snpsift/snpsift.yaml -p ${conda_dir}/conda_snpsift/snpsift_v5.1 ; fi


## Annotate vcf file using snpEff -noInteraction -noNextProt

mkdir -p ${out_dir}

# snpEFF

source activate ${conda_dir}/conda_snpeff/snpeff_v5.1

snpEff -Xmx10g ${assembly} -noInteraction -noNextProt -nodownload -dataDir ${snpeff_db} ${vcf} | gzip > ${out_dir}/${name}_snpEff.vcf.gz

source deactivate

# SnpSift

## Annotate using known databases such as Cosmic, ICGC (cbioportail) & cancerHotspots & gnomad

source activate ${conda_dir}/conda_snpsift/snpsift_v5.1

# Cosmic 

ANN_DB="${in_dir}/${assembly}/CosmicCodingMuts_v95_${assembly}.vcf.gz"
TAG_DB="COSMIC"

SnpSift -Xmx10g annotate -tabix -noId -noInfo -exists ${TAG_DB} ${ANN_DB} ${out_dir}/${name}_snpEff.vcf.gz | gzip > ${out_dir}/${name}_snpEff.Cosmic.vcf.gz

# ICGC

# Doesn't work if the database is gzipped 
ANN_DB="${in_dir}/${assembly}/icgc_release28_summary_simple_somatic_mutation.aggregated_${assembly}.vcf"
TAG_DB="ICGC"

SnpSift -Xmx10g annotate -noId -noInfo -exists ${TAG_DB} ${ANN_DB} ${out_dir}/${name}_snpEff.Cosmic.vcf.gz | gzip > ${out_dir}/${name}_snpEff.Cosmic.ICGC.vcf.gz

# CancerHotspot

ANN_DB="${in_dir}/${assembly}/cancerhotspots.v2.${assembly}.sort.vcf.gz"
TAG_DB="CancerHotspot"

SnpSift -Xmx10g annotate -tabix -noId -noInfo -exists ${TAG_DB} ${ANN_DB} ${out_dir}/${name}_snpEff.Cosmic.ICGC.vcf.gz | gzip > ${out_dir}/${name}_snpEff.Cosmic.ICGC.CancerHotspot.vcf.gz

# GnomAD

#gunzip -c ${in_dir}/${assembly}/af-only-gnomad.${assembly}.vcf.gz > ${in_dir}/${assembly}/af-only-gnomad.${assembly}.vcf

# Doesn't work if the database is gzipped 
ANN_DB="${in_dir}/${assembly}/af-only-gnomad.${assembly}.vcf"
TAG_DB="AF"

SnpSift -Xmx10g annotate -noId -Info ${TAG_DB} ${ANN_DB} ${out_dir}/${name}_snpEff.Cosmic.ICGC.CancerHotspot.vcf.gz | gzip > ${out_dir}/${name}_snpEff.Cosmic.ICGC.CancerHotspot.gnomad.vcf.gz

# dbNSFP

# Annotate using dbNSFP for SIFT/Polyphen/CADD 
ANN_DB="${in_dir}/${assembly}/dbNSFP4.1a_${assembly}.txt.gz"
FDB="SIFT_pred,Polyphen2_HDIV_pred,CADD_phred"

SnpSift dbnsfp -db ${ANN_DB} -f ${FDB} -collapse ${out_dir}/${name}_snpEff.Cosmic.ICGC.CancerHotspot.gnomad.vcf.gz | gzip > ${out_dir}/${name}_snpEff.Cosmic.ICGC.CancerHotspot.gnomad.dbNSFP.vcf.gz

# Split multiple effect in annotation into several lines

gunzip -c ${out_dir}/${name}_snpEff.Cosmic.ICGC.CancerHotspot.gnomad.dbNSFP.vcf.gz | ${conda_dir}/conda_snpsift/snpsift_v5.1/share/snpsift-5.1-0/scripts/vcfEffOnePerLine.pl | gzip > ${out_dir}/${name}_snpEff.Cosmic.ICGC.CancerHotspot.gnomad.dbNSFP.onePerLine.vcf.gz


# Extract defined fields into tabulated file

dbVar=$(for i in $(echo ${FDB} | sed "s|,| |g"); do echo "dbNSFP_"${i}; done)

SAMPLES=$(gunzip -c ${out_dir}/${name}_snpEff.Cosmic.ICGC.CancerHotspot.gnomad.dbNSFP.onePerLine.vcf.gz | grep CHROM | sed -n 's/^.*FORMAT//p')

var="" ; for j in ${SAMPLES} ; do var=${var}" "$(echo -e "GEN[${j}].AD"); done

FIELDS="CHROM POS REF ALT ANN[0].GENE ANN[0].IMPACT ANN[0].EFFECT ANN[0].FEATUREID ANN[0].RANK ANN[0].HGVS_C ANN[0].HGVS_P AF COSMIC ICGC CancerHotspot ${dbVar} ${var}"

SnpSift extractFields -e "."  -s ";" ${out_dir}/${name}_snpEff.Cosmic.ICGC.CancerHotspot.gnomad.dbNSFP.onePerLine.vcf.gz ${FIELDS}  > ${out_dir}/${name}_snpEff.Cosmic.ICGC.CancerHotspot.gnomad.dbNSFP.tsv

# Discard non coding variant "NR_"

grep -v "NR_" ${out_dir}/${name}_snpEff.Cosmic.ICGC.CancerHotspot.gnomad.dbNSFP.tsv | sed "s|ANN\[0\].||g" | sed "s|GEN\[||g" | sed "s|\]||g" | sed "s|dbNSFP_||g"  > ${out_dir}/${name}_snpEff.Cosmic.ICGC.CancerHotspot.gnomad.dbNSFP.coding.tsv

source deactivate



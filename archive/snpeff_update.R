# Construct file containing genes available in oncokb, mane and canonical files
# /bioinfo/local/build/Centos/R/R-4.1.0/bin/R

if (!require("tidyverse")) {install.packages("tidyverse", dependencies=TRUE, repos='http://cloud.r-project.org/')}

suppressMessages(require(tidyverse))

# no filtering 
min_freq <- 0
min_dp <- 0 

args <- commandArgs(trailingOnly = TRUE) 

in_dir <- args[1] ; stopifnot(!is.na(in_dir)) 
tsv_file <- args[2] ; stopifnot(!is.na(tsv_file)) 
tumor_name <- args[3] ; stopifnot(!is.na(tumor_name))  
assembly <- args[4] ; stopifnot(!is.na(assembly))  
out_dir <- args[5] ; stopifnot(!is.na(out_dir))  
min_freq <- as.numeric(args[6])
min_dp <- as.numeric(args[7])


dir.create(out_dir,showWarnings=FALSE)

# in_dir <- "/data/kdi_prod/.kdi/project_workspace_0/1198/acl/02.00/snp_eff/databases"
# tsv_file <- paste0(in_dir,"/snpeff/results/D777E181_vs_D777E133_snpEff.Cosmic.ICGC.CancerHotspot.gnomad.dbNSFP.coding.tsv")
# tumor_name <- "D777E181"
# assembly <- "hg38" 
# out_dir <- "/data/kdi_prod/.kdi/project_workspace_0/1198/acl/02.00/snp_eff/databases/snpeff/results"

# min_freq <- 5 # set at 0 to remove filtering
# min_dp <- 20 # set at 0 to remove filtering


# Load  NM with longest size and shortest name

hg38_nm_length_file <- paste0(in_dir,"/hg38/hg38_gencode_v39_genepred.length_nm.txt")
hg19_nm_length_file <- paste0(in_dir,"/hg19/hg19_gencode_v39lift37_genepred.length_nm.txt")

hg38_nm_length <- read.table(hg38_nm_length_file, header = TRUE, fill = TRUE)
hg19_nm_length <- read.table(hg19_nm_length_file, header = TRUE, fill = TRUE)

# Load transcripts of interest (oncokb,mane, canonical)

transcripts_int_file <- paste0(in_dir,"/oncokb_mane_canonical_transcript.txt")
t_nm <- read.table(transcripts_int_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)

## Select variants of interest

variant_file <- paste0(in_dir,"/snpeff_genetic_variant_effect.csv")
variant <- read.table(variant_file,sep=",",header=TRUE, check.names=FALSE)

# Load annotated tsv file

tsv_onco_oncoKB_file <- paste0(out_dir,"/",tumor_name,"_oncoprint_oncoKB_table.tsv")
tsv_onco_oncoKB_merge_file <- paste0(out_dir,"/",tumor_name,"_oncoprint_oncoKB_merge_table.tsv")
tsv_onco_merge_file <- paste0(out_dir,"/",tumor_name,"_oncoprint_merge_table.tsv")
tsv_onco_file <- paste0(out_dir,"/",tumor_name,"_oncoprint_table.tsv")
tsv_partial_file <- paste0(out_dir,"/",tumor_name,"_partial_table.tsv")
tsv_full_file <- paste0(out_dir,"/",tumor_name,"_full_table.tsv")

print(tsv_partial_file)


#### 

tsv <- read.table(tsv_file, header=TRUE, check.names=FALSE)

tsv1 <- tsv %>% select(chr = CHROM, pos = POS, ref = REF, alt = ALT, Gene = GENE, Impact = IMPACT, Effect = EFFECT, Transcript = FEATUREID, rank = RANK, "c." = HGVS_C, "p." = HGVS_P, AF, COSMIC, ICGC, CancerHotspot, SIFT = SIFT_pred,Polyphen_HDIV = Polyphen2_HDIV_pred, CADD_phred, ends_with(".GT"), ends_with(".AD"))

tsv2 <- tsv1 %>% mutate(First_Effect=str_remove(Effect,"&.*$"), 
                        NM=str_remove(Transcript,"\\..*$"),
                        Databases=ifelse(COSMIC == "true"| ICGC == "true" | CancerHotspot == "True", "Yes", "No")) %>% 
                left_join(t_nm[,grep("is|NM_onco",colnames(t_nm))], by=c("NM"="NM_onco"))

tsv3 <- tsv2 %>% left_join(variant, by=c("First_Effect"="snpEff")) %>% 
                 mutate(Decision=case_when(Genetic == "TERT" & Gene == "TERT" & (str_detect(c., "c.-(146|124|57|138|139)[ATCG]") | str_detect(c.,"-139(ins|del)")) ~ "TERT",
                                       Genetic == "TERT" ~ NA_character_,  
                    TRUE ~ Genetic))

# just C>T, also in -57T>G , might also be a double mutation next to each other but often called in two different lines
# -(124|146)((ins|del)|[ACTG])
 

tsv4 <- tsv3  %>% separate(glue::glue("{tumor_name}.AD"), into=c("RD","AD"),sep=",") %>% mutate(!!glue::glue("{tumor_name}.DP") :=  as.numeric(RD)+as.numeric(AD)) %>% mutate(!!glue::glue("{tumor_name}.FREQ") :=  round(as.numeric(AD)*100/(as.numeric(RD)+as.numeric(AD)),2)) 


# add informations about polyphen & sift 

tsv5 <- tsv4  %>% mutate(SIFT_bis= gsub(".,","",gsub(",.","",SIFT,fixed=TRUE),fixed=TRUE), Polyphen_bis=gsub(".,","",gsub(",.","",Polyphen_HDIV,fixed=TRUE),fixed=TRUE)) %>% mutate(SIFT_prediction = case_when(SIFT_bis=="D" ~ "Deleterious", SIFT_bis=="T" ~ "Tolerated", SIFT_bis %in% c("D,T", "T,D") ~ "Ambiguous", TRUE ~ NA_character_), Polyphen_prediction = case_when(Polyphen_bis=="D" ~ "Deleterious", Polyphen_bis=="B" ~ "Benign", Polyphen_bis=="P" ~ "Probably damaging", Polyphen_bis %in% c("P,D","D,P","B,D","D,B","B,P","P,B","D,B,P", "D,P,B", "B,D,P", "B,P,D", "P,D,B","P,B,D") ~ "Ambiguous", TRUE ~ NA_character_))  %>% select(-SIFT_bis,-Polyphen_bis) 

write.table(tsv5 %>% dplyr::select(-c(Oncogene:Signature)) %>% mutate(t=as.numeric(str_remove(str_replace(str_replace(chr,"chrX","chr23"), "chrY","chr24"), "chr"))) %>% arrange(t, pos) %>% select(-t), tsv_full_file, col.names=TRUE, row.names=FALSE, sep="\t",quote=FALSE)

# # find nm length 

# assign(glue::glue("{assembly}_nm_length_max"),get(glue::glue("{assembly}_nm_length")) %>% filter(grepl("NM",refseq_mrna)) %>% filter(refseq_mrna %in% tsv5$NM) %>% group_by(hgnc_symbol) %>% slice_max(glue::glue("length_{assembly}")) %>%  ungroup())

# hard filtering based on dp and freq of the tumor & filter on genetic effect (to check)
tsv6 <- tsv5 %>% filter(.data[[glue::glue("{tumor_name}.DP")]]>=min_dp & .data[[glue::glue("{tumor_name}.FREQ")]]>=min_freq & !is.na(Decision)) 

tsv6b <- tsv6 %>% left_join((get(glue::glue("{assembly}_nm_length"))[,c("refseq_mrna",glue::glue("length_{assembly}"))]  %>% distinct()), by=c("NM"="refseq_mrna"))

# select nm in oncokb > mane then > canonical then longest transcript then (if several), the one with the shortest NM (found in an early stage)

tsv7 <- tsv6b %>% mutate(Selection=case_when(NM %in% t_nm$NM_onco ~ "ONCOKB", NM %in% t_nm$NM_Mane ~ "Mane", NM %in% t_nm[,glue::glue("NM_{assembly}_canonical")] ~ "Canonical" ,TRUE ~ NA_character_), Selection = factor(Selection, levels = c("ONCOKB", "Mane", "Canonical")))

tsv8 <- tsv7 %>%  mutate(var = paste0(chr,pos,ref,alt), NM_number=as.numeric(str_remove(NM,"NM_"))) %>% group_by(var) %>% arrange(Selection, -length_hg38, NM_number) %>% slice_head() %>% ungroup() %>% select(-var) %>% mutate(Selection = replace_na(as.character(Selection), "Longest_size"))   %>% as.data.frame() 

write.table(tsv8 %>% dplyr::select(-c(Oncogene:Signature, length_hg38, Selection, NM_number)) %>% mutate(t=as.numeric(str_remove(str_replace(str_replace(chr,"chrX","chr23"), "chrY","chr24"), "chr"))) %>% arrange(t, pos) %>% select(-t), tsv_partial_file, col.names=TRUE, row.names=FALSE, sep="\t",quote=FALSE)

## table report for oncoprint : 
# keep variant with databases & # Keep variant with yes for oncoprint & 

tsv_onco <- tsv8 %>% filter(Oncoprint %in% c("Yes", "TERT"))  %>% 
mutate(Onco_status = ifelse(
((grepl("databases",Oncogene) & Databases == "Yes" & isOncogene == "Yes") | (grepl("databases",TSG) & Databases == "Yes" & isTSG == "Yes") | (grepl("databases",Both) & Databases == "Yes" & isBoth == "Yes")) | 
(Oncoprint == "Yes" & is.na(isOncoKB)) |
((Oncogene =="TERT"   & isOncogene == "Yes") | (TSG == "Yes" & isTSG == "Yes") | (Both == "Yes" & isBoth == "Yes")), "Oncoprint",NA_character_
))

tsv_onco <- tsv_onco %>% filter(!is.na(Onco_status))
write.table(tsv_onco %>%  dplyr::select(-c(Oncogene:Signature, length_hg38, Selection, NM_number,  Onco_status)) %>% mutate(t=as.numeric(str_remove(str_replace(str_replace(chr,"chrX","chr23"), "chrY","chr24"), "chr"))) %>% arrange(t, pos) %>% select(-t), tsv_onco_file, col.names=TRUE, row.names=FALSE, sep="\t",quote=FALSE)

tsv_onco_merge <- tsv_onco %>% group_by(Gene) %>% arrange(Decision) %>% summarize(!!glue::glue("{tumor_name}") := paste(unique(Decision), collapse=";")) %>%  ungroup() 
write.table(tsv_onco_merge, tsv_onco_merge_file, col.names=TRUE, row.names=FALSE, sep="\t",quote=FALSE)


tsv_onco_oncoKB <- tsv_onco %>% filter(!is.na(isOncoKB)) 
write.table(tsv_onco_oncoKB %>%  dplyr::select(-c(Oncogene:Signature, length_hg38, Selection, NM_number, Onco_status)) %>% mutate(t=as.numeric(str_remove(str_replace(str_replace(chr,"chrX","chr23"), "chrY","chr24"), "chr"))) %>% arrange(t, pos) %>% select(-t), tsv_onco_oncoKB_file, col.names=TRUE, row.names=FALSE, sep="\t",quote=FALSE)

tsv_onco_oncoKB_merge <- tsv_onco_oncoKB %>% group_by(Gene) %>% arrange(Decision) %>% summarize(!!glue::glue("{tumor_name}") := paste(unique(Decision), collapse=";")) %>%  ungroup() 
write.table(tsv_onco_oncoKB_merge, tsv_onco_oncoKB_merge_file, col.names=TRUE, row.names=FALSE, sep="\t",quote=FALSE)

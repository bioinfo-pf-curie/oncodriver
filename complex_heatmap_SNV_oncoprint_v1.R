#SNV
in.dir <- ""
files <- dir(path = in.dir, pattern = "*_oncoprint_oncoKB_table.tsv")
snv.data <- files %>%
  map_dfr(~ {
    tumor <- gsub("_oncoprint_oncoKB_table.tsv","", basename(.))
    read.table(file.path(in.dir, .), header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE) %>% 
      mutate(across(.fns = as.character)) %>% mutate(pos = as.numeric(pos)) %>% 
      mutate(ID = tumor) %>% 
      select(-contains(c(".AD", ".DP", ".FREQ", "CADD_phred"))) %>% 
      mutate(AF = na_if(AF, ".")) %>% mutate(AF = as.numeric(AF)) %>%  filter(AF < 0.001 | is.na(AF)) %>%
      mutate(end = pos + 1) %>% 
      rename(chrom = chr, start  = pos) %>%  relocate(chrom, start, end)
  }) 

# loh.intersect <- snv.data %>% group_by(ID) %>% bed_intersect(group_by(all.cnv.files, ID)) %>% 
#   select(ID.x, ID.y, LOH.y, chrom, start.x, end.x) 

mat <- snv.data %>% 
  #full_join(., y=all.cnv.oncokb, ID, by = c("ID", "Gene", "chrom", "start", "end", "isOncogene" = "Is.Onco", "isTSG" = "Is.TSG", "isBoth" = "Is.Both", "Decision" = "call")) %>% 
  #full_join(., y=loh.intersect, by = c("ID"="ID.x", "chrom", "start"="start.x", "end"="end.x")) %>% #select(-c(LOH, ID.y)) %>% rename(LOH = LOH.y) %>% 
  select(ID, chrom, start, end, Gene, Decision, isOncogene, isTSG, isBoth, isUnknown, COSMIC, ICGC, CancerHotspot, LOH) %>% 
  filter( 
    ((isOncogene == "Yes" & isTSG == "No") & (
      (Decision == "AMP") |
        (Decision == "inframe indel" & (COSMIC == "true" | ICGC == "true" | CancerHotspot == "true") ) |
        (Decision == "missense" & (COSMIC == "true" | ICGC == "true" | CancerHotspot == "true") ))) |
      
      (Decision == "TERT") |
      
      ( (isTSG == "Yes" & isOncogene == "No") & (
        (Decision == "DEL") | (Gene == "POLE") | (Gene == "POLD1") |
          (Decision == "inframe indel" & (COSMIC == "true" | ICGC == "true" | CancerHotspot == "true")) |
          (Decision == "stop gain") |
          (Decision == "frameshift indel") |
          (Decision == "splicing") |
          (Decision == "missense" & (COSMIC == "true" | ICGC == "true" | CancerHotspot == "true")) #& LOH == 1) 
      )) |
      
      ( (isBoth == "Yes" | (isOncogene == "No" & isTSG == "No")) & (
        (Decision == "AMP") |
          (Decision == "DEL") |
          (Decision == "inframe indel" & (COSMIC == "true" | ICGC == "true" | CancerHotspot == "true")) |
          (Decision == "stop gain") | 
          (Decision == "frameshift indel") | 
          (Decision == "splicing") |
          (Decision == "missense" & (COSMIC == "true" | ICGC == "true" | CancerHotspot == "true")) ))
  ) %>% 
  select(Gene, ID, Decision) %>%
  group_by(Gene, ID) %>%
  distinct() %>% 
  mutate(Decision = paste0(Decision, collapse = ";")) %>%
  distinct() %>%
  ungroup() %>% 
  pivot_wider(.,names_from = ID, values_from = Decision) %>% 
  column_to_rownames(var = "Gene") %>% 
  as.data.frame()

# filtered.samples <- patientID[which(!patientID$SampleID %in% colnames(mat)), "SampleID"]
# mat[,filtered.samples] <- NA

mat <- as.matrix(mat)
mat[is.na(mat)] <- ""


### Oncoprint
fq <- apply(mat, 1, function (x) sum((x != ""), na.rm = TRUE) / ncol(mat))
fqOrdered <- fq[order(fq, decreasing = TRUE)]
newmat <- mat[which(fq >= 0.05),]
new_order <- names(fqOrdered[(fqOrdered >= 0.05)])

col = c("background"="#CCCCCC",
        "frameshift indel" ="orange", 
        "inframe indel"="yellowgreen", 
        "splicing"="hotpink", #"pink", 
        unknown_splice = "darkorchid4",
        "missense"="forestgreen", 
        synonymous_variant="mediumorchid2",
        "stop gain"="gold", 
        start_lost = "yellow2", 
        stop_lost = "red",
        "unknown"="darkgrey",
        DEL = 'deepskyblue',
        'AMP' = 'red3',
        'TERT' = 'dodgerblue4')

################################################################################
### Oncoprint OncoKB
################################################################################

# oncogene <- read.csv("oncoKB_cancerGeneList.tsv", header = TRUE, sep = "\t")
# oncogene <- mutate(oncogene, class = case_when(
#   Is.Oncogene == "Yes" & Is.Tumor.Suppressor.Gene == "Yes"  ~ "Both",
#   Is.Oncogene == "Yes" & Is.Tumor.Suppressor.Gene == "No"  ~ "Oncogene",
#   Is.Oncogene == "No" & Is.Tumor.Suppressor.Gene == "Yes"  ~ "TSG",
#   Is.Oncogene == "No" & Is.Tumor.Suppressor.Gene == "No"  ~ "Unknown",
# ))
# 
# oncogene$class <- factor(oncogene$class, levels = c("Oncogene", "TSG", "Both", "Unknown"))
# oncogenemat <- newmat[which(rownames(newmat) %in% oncogene$Hugo.Symbol),]
# oncogene_order <- names(fqOrdered[which(names(fqOrdered) %in% rownames(oncogenemat))])

heatmap_legend_param = list(title = "Alterations", at = c("AMP", "DEL", "TERT", "stop gain", "frameshift indel", "missense", "splicing"), 
                            labels = c("AMP", "DEL", "Upstream gene variant (only TERT)", "Stop gain", "Frameshift indel", "Missense", "Splicing"),
                            ncol = 1,
                            legend_gp = gpar(fontsize = 8))

op1 <- oncoPrint(newmat,
                 alter_fun = list(
                   background = function(x, y, w, h) 
                     grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#CCCCCC", col = NA)),
                   
                   
                   AMP = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, 
                                                        gp = gpar(fill = col["AMP"], col = NA)), 
                   
                   DEL = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, 
                                                        gp = gpar(fill = col["DEL"], col = NA)),
                   
                   TERT = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.8, 
                                                         gp = gpar(fill = col["TERT"], col = NA)), 
                   
                   `stop gain` = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.8, 
                                                                gp = gpar(fill = col["stop gain"], col = NA)),
                   
                   `frameshift indel` = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.65, 
                                                                       gp = gpar(fill = col["frameshift indel"], col = NA)),
                   
                   `inframe indel` = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.6, 
                                                                    gp = gpar(fill = col["inframe indel"], col = NA)),
                   
                   missense = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.5, 
                                                             gp = gpar(fill = col["missense"], col = NA)),
                   
                   splicing = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.3, 
                                                             gp = gpar(fill = col["splicing"], col = NA)),
                   
                   stop_lost = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.2, 
                                                              gp = gpar(fill = col["stop_lost"], col = NA)),
                   
                   start_lost = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.2, 
                                                               gp = gpar(fill = col["start_lost"], col = NA))
                   
                   
                   
                 ), col = col, 
                 #top_annotation = column_ha,
                 # row_order = oncogene_order,
                 row_names_gp = gpar(fontsize = 8), pct_gp = gpar(fontsize = 5),
                 # row_split = oncogene[match(rownames(oncogenemat), oncogene$Hugo.Symbol),"class"],
                 gap = unit(c(5),"mm"),
                 remove_empty_columns = FALSE, show_column_names = FALSE, alter_fun_is_vectorized = TRUE) 


# jpeg("Oncoprint_CNV_SNV_OncoKB_genes_freq3_23-06-05.jpeg",
#      units="in", width=12, height=7, res = 600)
draw(op1, heatmap_legend_side = "bottom")#, annotation_legend_side = "bottom")
# dev.off()

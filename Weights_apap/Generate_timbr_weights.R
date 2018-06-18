# Pannala et al., (2018)
# ################################################################################### #
# ############################### Load packages  #################################### #
# ################################################################################### #
# inatall appropriate packages in R to load the librarys
options(java.parameters = "-Xmx4g")
options(stringsAsFactors = FALSE)
source("ncomm_helper.R")

library(limma)
library(Biobase)
library(reshape2)
library(dplyr)
library(ggplot2)
library(xlsx)
# provide path to the current directory
path.manuscript.directory = "(Insert_current_directory)/"
path.expression.directory = "(Provide path location of the sleuth tables)/"
path.weights.directory = "(Provide path location to store the generated weights)/"

rxn.info.load = read.xlsx2(paste0(path.manuscript.directory,"Supplementary_iRno_v2.xlsx"),
                           sheetIndex = 1, startRow = 2) %>% as.tbl

rxn.info = rxn.info.load %>% 
  mutate(enabled = ifelse(!is.na(enabled), enabled == "true", F)) %>% 
  mutate(enabled_hsa = ifelse(!is.na(enabled_hsa), enabled_hsa == "true", F)) %>% 
  mutate(enabled_rno = ifelse(!is.na(enabled_rno), enabled_rno == "true", F))%>% 
  mutate(n_hsa = as.numeric(as.character(n_gene_hsa)),
         n_rno = as.numeric(as.character(n_gene_rno)))

rxn.gene = rxn.info %>% filter(enabled) %>% 
  select(rxn_id, hsa = gpr_hsa, rno = gpr_rno) %>%
  melt(c("rxn_id")) %>% ef_df %>% 
  # Excel annoyingly transforms GPR rules into date/time values
  mutate(value = gsub("193.5625","4645:30",value,fixed = T)) %>%
  mutate(value = gsub("[\\(\\)\\;\\:]+",";",value)) %>%
  ef_split_value(";") %>% ef_df %>% 
  rename(organism_id = variable, gene_id = value) %>% 
  filter(nchar(gene_id) > 0) %>% distinct 

# all gene_id values should be integers
rxn.gene %>% count(grepl("^[0-9]+$",gene_id))

# Conversion tables for rat genomes from Ensembl/Biomart
# Used to convert Ensembl transcript IDs (target_id) to NCBI/Entrez gene IDs (gene_id)
# Biomart_rno_ids.txt.gz can be downloaded from the Ensemble biomart by choosing appropriate annotations
# for the latest release Ensemble_gene_version and select Rnor_version using the link below
# http://useast.ensembl.org/biomart/martview/
biomart.ncbi.rno.load = read.table("Biomart_rno_ids.txt.gz", 
                                   sep = "\t", quote = "", header = T) %>% as.tbl %>%
  setNames(c("ensembl_gene_id", "ensembl_transcript_id", "ensembl_gene_type", "ensembl_transcript_type",
             "ensembl_gene_status", "ensembl_transcript_status", "ensembl_chromosome", "ensembl_description",
             "ensembl_gene_name", "ensembl_transcript_name","gene_id")) %>% 
  mutate(gene_id = ifelse(!is.na(gene_id), as.character(gene_id), "")) %>%
  mutate(organism_id = "rno")

# Load expression changes from sleuth (alterantively one can change this for 10 h results)
differential.expression.file = paste0(path.expression.directory, "Sleuth_table_liver_5h.csv")
differential.expression.load = read.csv(differential.expression.file, header = T, stringsAsFactors = F) %>% as.tbl
differential.expression.data = differential.expression.load %>%
  mutate(organism_id = "rno", drug_id = "acetaminophen", dose_id = "d1", time_id = "t1") %>%
  mutate(efit_id = paste0(organism_id,"_",drug_id,"_",time_id,"_",dose_id)) %>%
  mutate(efit_root = paste0("dod","_","rnaseq","_", organism_id,"_",drug_id,"_",time_id,"_",dose_id)) %>%
  mutate(target_id = as.character(target_id), feature_id = target_id, ensembl_transcript_id = target_id) %>% 
  mutate(fdr = qval, logfc = log2(exp(b)), ave = mean_obs) %>%
  filter(!is.na(logfc), !is.na(fdr)) %>%
  inner_join(biomart.ncbi.rno.load, by = c("organism_id","ensembl_transcript_id")) %>%
  arrange(desc(var_obs), fdr, pval, desc(abs(logfc))) %>% 
  group_by(efit_id, efit_root, organism_id, drug_id, time_id, dose_id, gene_id) %>% slice(1) %>% ungroup

# check unique rat genes map to this RNA-seq dataset with fdr/logfc values
differential.expression.data %>% 
  select(organism_id, gene_id) %>% distinct %>%
  count(organism_id, gene_id %in% c(rxn.gene[["gene_id"]]))

fdr.cutoff = 0.1

# Check unique rat genes were significantly differentially expressed
differential.expression.data %>% filter(fdr < fdr.cutoff) %>%
  select(organism_id, gene_id) %>% distinct %>%
  count(organism_id, gene_id %in% c(rxn.gene[["gene_id"]]))


metabolic.differential.expression = differential.expression.data %>%
  semi_join(rxn.gene %>% select(gene_id, organism_id) %>% distinct) %>% 
  mutate(significant = fdr < fdr.cutoff,
         direction = ifelse(significant, sign(logfc),0)) %>%
  group_by(efit_id, efit_root, organism_id, drug_id, time_id, dose_id) %>%
  mutate(model_data_count = n(),
         model_feature_count = length(unique(feature_id)),
         model_gene_count = length(unique(gene_id)),
         n_up = sum(direction > 0), 
         n_dn = sum(direction < 0), 
         n_significant = sum(direction != 0)) %>% ungroup %>%
  mutate(pct_significant = 100 * n_significant / model_gene_count) %>%
  mutate(efit_ok = (model_gene_count == model_data_count) & (pct_significant > 1))

rxn.pubmed = rxn.info %>% filter(enabled) %>% select(rxn_id, pubmed_id) %>% 
  melt(c("rxn_id")) %>% ef_df %>% 
  mutate(value = gsub("\\-","",value)) %>% 
  mutate(value = gsub("PMID[\\:\\-]*",";PMID:",value)) %>% 
  ef_split_value(";") %>% ef_df %>% 
  filter(grepl("PMID|DOI|UNIPROT", value)) %>%
  filter(grepl("[0-9]+",value)) %>% distinct

rxn.gpr.size = bind_rows(list(
  rxn.gene %>%
    mutate(variable = paste0(organism_id, "_count_all")) %>% 
    group_by(rxn_id, variable) %>% 
    summarize(value = length(unique(setdiff(gene_id,"0")))) %>% ungroup,
  rxn.gene %>%
    # semi_join(metabolic.differential.expression %>% select(organism_id, gene_id) %>% distinct) %>% 
    left_join(metabolic.differential.expression %>% 
                select(organism_id, gene_id) %>% distinct %>% mutate(gene_overlap = T)) %>% 
    mutate(gene_overlap = ifelse(!is.na(gene_overlap) & !grepl("^0$", gene_id), gene_overlap, F)) %>%
    mutate(variable = paste0(organism_id, "_count")) %>% 
    group_by(rxn_id, variable) %>% 
    summarize(value = length(unique(setdiff(gene_id,"0")))) %>% ungroup)) %>%
  dcast(rxn_id ~ variable, value.var = "value", fill = 0) %>% ef_df

timbr.weights.default = rxn.info %>% filter(enabled) %>% 
  select(rxn_id,rxn_class) %>% distinct %>% 
  left_join(rxn.pubmed %>% count(rxn_id) %>% ungroup %>% rename(pubmed_count = n)) %>% 
  mutate(pubmed_count = ifelse(!is.na(pubmed_count), pubmed_count, 0)) %>% 
  left_join(rxn.gpr.size) %>% 
  mutate(rno_count = ifelse(!is.na(rno_count), rno_count, 0),
         rno_count_all = ifelse(!is.na(rno_count_all), rno_count_all, 0),
         hsa_count = ifelse(!is.na(hsa_count), hsa_count, 0),
         hsa_count_all = ifelse(!is.na(hsa_count_all), hsa_count_all, 0)) %>%
  mutate(rxn_enzymatic = hsa_count > 0 | rno_count > 0,
         rxn_referenced = pubmed_count > 0 ) %>%
  mutate(timbr_weight_enzymatic = ifelse(rxn_enzymatic,1,2),
         timbr_weight_referenced = ifelse(rxn_referenced,1,2),
         timbr_weight_class = ifelse(rxn_class == "boundary",2, ifelse(rxn_class == "transport",2,1)),
         timbr_weight_default = timbr_weight_enzymatic * timbr_weight_referenced * timbr_weight_class)

timbr.rxn.setup = timbr.weights.default %>% select(rxn_id, timbr_weight_default) %>%
  left_join(rxn.info %>% select(rxn_id, rno = gpr_rno, hsa = gpr_hsa)) %>%
  melt(c("rxn_id","timbr_weight_default")) %>% ef_df %>%
  mutate(value = gsub("193.5625","4645:30",value, fixed = T))

timbr.expression.setup = metabolic.differential.expression %>% 
  mutate(limma_id = paste0(organism_id, "_", drug_id, "_", time_id, "_", dose_id),
         limma_ok = efit_ok) %>% filter(limma_ok) %>% 
  ef_df_slice("limma_id")

timbr.weights.list = timbr.expression.setup %>% 
  lapply(ef_timbr_weights,timbr.rxn.setup, gene.ignore = c(), gene.keep.all = F)
timbr.weights = timbr.weights.list %>% bind_rows

timbr.weights %>% with(qplot(timbr_weight_ctl, timbr_weight_trt))

# Simulating TIMBR predictions requires an irreversible metabolic network, which 
# appends *_f or *_r to each reaction identifier in the forward or reverse direction
rxn.irreversible = bind_rows(list(
  timbr.weights %>% select(rxn_id) %>% distinct %>% mutate(irxn_id = paste0(rxn_id,"_f")),
  timbr.weights %>% select(rxn_id) %>% distinct %>% mutate(irxn_id = paste0(rxn_id,"_r"))))

timbr.weights.irreversible = bind_rows(list(
  timbr.weights %>% mutate(rxn_irreversible = paste0(rxn_id,"_f")),
  timbr.weights %>% mutate(rxn_irreversible = paste0(rxn_id,"_r"))))

# Organize reaction weights for each individual condition:
# (treatment/control, organism, compound, treatment dose, treatment duration)
rno.timbr.weights = bind_rows(list(
  timbr.weights.irreversible  %>%
    filter(organism_id == "rno") %>% 
    mutate(timbr_id = paste0(limma_id, "_ctl")) %>%
    select(timbr_id, organism_id, rxn_id, rxn_irreversible, rxn_weight = timbr_weight_ctl),
  timbr.weights.irreversible  %>%
    filter(organism_id == "rno") %>% 
    mutate(timbr_id = paste0(limma_id, "_trt")) %>%
    select(timbr_id, organism_id, rxn_id, rxn_irreversible, rxn_weight = timbr_weight_trt))) %>%
  dcast(organism_id + rxn_id + rxn_irreversible ~ timbr_id, value.var = "rxn_weight") %>% ef_df

# Save outputs in a folder that can be accessed by MATLAB
write.table(rno.timbr.weights, file = paste0(path.weights.directory,"timbr_weights_rno_liver5h.txt"), 
            sep = "\t", quote = F, row.names = F)

# Same script can be used for "Sleuth_table_liver_10h.csv" by changing the file names accordingly

rm(list = ls()) # R version 4.3.1 (2023-06-16)
library(tidyverse) # tidyverse_2.0.0
setwd(".")
out.dir <- "2.0/"

# --- Data ---
# PanDrugs table
pandrugs <- 
  read.table(list.files(out.dir, pattern = "PanDrugs_.*_gs_dr.tsv", full.names = TRUE), 
             header = TRUE, sep = "\t", quote = "")

# Salt synonyms
synonyms <- read.table("data/synonyms.tsv", header = TRUE, sep = "\t")

# --- Code ---
# Select a reference synonym for standard_drug_name
synonyms <- synonyms %>%
  mutate(reference = new_show_drug_name == show_drug_name) %>%
  group_by(new_show_drug_name) %>%
  mutate(has.reference = any(reference),
         label = seq(n())) %>%
  ungroup() %>%
  mutate(reference = if_else(!has.reference & label == 1, TRUE, reference)) %>%
  select(-c(has.reference, label))

# New PanDrugs table
pandrugs.new <- pandrugs

# Standardise names
pandrugs.new <- pandrugs.new %>%
  mutate(is.indatraline = show_drug_name == "(R,S)-INDATRALINE",
         is.conjugatedestrogen = source_drug_name == "CONJUGATED ESTROGENS",
         standard_drug_name = 
           case_when(is.indatraline ~ "INDATRALINE",
                     is.conjugatedestrogen ~ "CONJUGATED ESTROGENS",
                     TRUE ~ standard_drug_name),
         show_drug_name = 
           case_when(is.indatraline ~ "INDATRALINE",
                     is.conjugatedestrogen ~ "CONJUGATED ESTROGENS",
                     TRUE ~ show_drug_name)) %>%
  select(-is.indatraline, -is.conjugatedestrogen)

# Merge synonyms
pandrugs.new <- pandrugs.new %>%
  left_join(synonyms, by = "show_drug_name") %>%
  rename(ID = new_show_drug_name) %>%
  unique()

# Collapse drug info
updated <- pandrugs.new %>%
  filter(!is.na(ID)) %>%
  select(-show_drug_name) %>%
  unique() %>%
  group_by(ID) %>%
  mutate(standard_drug_name = unique(standard_drug_name[reference]),
         family = if_else(family == "Other", "1", family),
         family = paste(sort(unique(str_trim(unlist(
           str_split(family, pattern = ", "))))), collapse = ", "),
         family = str_remove(family, pattern = "^1, "),
         family = if_else(family == "1", "Other", family),
         n.status = length(unique(status)),
         is.approved = any(status == "Approved"),
         is.clinical = any(status == "Clinical trials"),
         is.experimental = any(status == "Experimental"),
         pathology = unique(pathology),
         is.cancer = any(cancer == "cancer"),
         is.clinicalcancer = any(cancer == "clinical cancer"),
         cancer = if_else(is.cancer | is.clinicalcancer | cancer == "", "1", cancer),
         cancer = paste(sort(unique(str_trim(unlist(
           str_split(cancer, pattern = "\\|"))))), collapse = " | "),
         cancer = str_remove(cancer, pattern = "^1 \\| "),
         cancer = if_else(cancer == "1", "", cancer),
         extra = if_else(extra == "", "1", extra),
         extra = paste(sort(unique(str_trim(unlist(
           str_split(extra, pattern = ", "))))), collapse = ", "),
         extra = str_remove(extra, pattern = "^1, "),
         extra = if_else(extra == "1", "", extra),
         extra2 = str_remove(paste(unique(extra2), collapse = "|"), 
                             pattern = "\\|")) %>%
  ungroup() %>%
  mutate(status = case_when(n.status > 1 & is.approved ~ "Approved",
                            n.status > 1 & is.clinical ~ "Clinical trials",
                            n.status > 1 & is.experimental ~ "Experimental",
                            n.status == 1 ~ status),
         cancer = case_when(status == "Approved" & !is.clinicalcancer ~ cancer,
                            status == "Approved" & is.clinicalcancer ~ 
                              "clinical cancer",
                            status == "Clinical trials" & is.cancer ~ "cancer",
                            status == "Clinical trials" ~ "",
                            status == "Experimental" ~ "")) %>%
  select(-c(reference:is.clinicalcancer)) %>%
  unique()

# Collapse drug-gene info
updated <- updated %>%
  group_by(gene_symbol) %>%
  mutate(checked_gene_symbol = unique(checked_gene_symbol)) %>%
  ungroup() %>%
  group_by(ID, checked_gene_symbol) %>%
  mutate(pathways = unique(pathways),
         n.target_marker = length(unique(target_marker)),
         ind_pathway = unique(ind_pathway),
         gene_dependency = unique(gene_dependency),
         gscore = unique(gscore),
         driver_gene = unique(driver_gene)) %>%
  ungroup() %>%
  mutate(target_marker = if_else(n.target_marker == 2, "target",
                                 target_marker)) %>%
  unique() %>%
  group_by(ID, source, checked_gene_symbol, resistance) %>%
  mutate(alteration = unique(alteration)) %>%
  ungroup() %>%
  group_by(ID, source, checked_gene_symbol, alteration) %>%
  mutate(resistance = unique(resistance)) %>%
  ungroup() %>%
  select(-n.target_marker) %>% 
  unique()

# Merge all entries and recompute DScores
updated <- updated %>%
  mutate(is.cancer = case_when(cancer %in% c("", "clinical cancer") ~ cancer,
                               extra == "solid tumors" ~ "cancer",
                               TRUE ~ "cancer"),
         dscore = 
           case_when(target_marker == "target" & status == "Approved" & 
                       is.cancer == "cancer" ~ 1,
                     target_marker == "marker" & status == "Approved" &
                       is.cancer == "cancer" ~ 0.9,
                     target_marker == "target" & status == "Approved" &
                       is.cancer == "clinical cancer" ~ 0.8,
                     target_marker == "marker" & status == "Approved" &
                       is.cancer == "clinical cancer" ~ 0.7,
                     target_marker == "target" & status == "Clinical trials" &
                       is.cancer == "cancer" ~ 0.6,
                     target_marker == "marker" & status == "Clinical trials" &
                       is.cancer == "cancer" ~ 0.5,
                     target_marker == "target" & status == "Approved" &
                       is.cancer == "" ~ 0.4,
                     target_marker == "marker" & status == "Approved" &
                       is.cancer == "" ~ 0.3,
                     target_marker == "target" & status == "Clinical trials" &
                       is.cancer == "" ~ 0.2,
                     target_marker == "marker" & status == "Clinical trials" &
                       is.cancer == "" ~ 0.1,
                     target_marker == "target" & status == "Experimental" &
                       is.cancer == "" ~ 0.0008,
                     target_marker == "marker" & status == "Experimental" &
                       is.cancer == "" ~ 0.0004),
         dscore = if_else(resistance == "resistance", -1 * dscore, dscore)) %>%
  rename(show_drug_name = ID) %>%
  select(any_of(colnames(pandrugs.new))) %>%
  mutate(version = "New")

# Add to PanDrugs table
pandrugs.new <- pandrugs.new %>%
  filter(is.na(ID)) %>%
  select(-c(ID, reference)) %>%
  mutate(version = "Old") %>%
  bind_rows(updated) %>%
  unique()

# Clean VITAMIN E entry (n.status == 2 & n.cancer == 2)
vitamin.e <- pandrugs.new %>%
  filter(show_drug_name == "VITAMIN E")

vitamin.e %>%
  filter(version == "Old") %>%
  select(standard_drug_name, status, cancer, dscore) %>%
  unique()

vitamin.e <- vitamin.e %>%
  group_by(show_drug_name) %>%
  mutate(standard_drug_name = "VITAMIN E",
         family = unique(family),
         status = "Withdrawn",
         pathology = unique(pathology),
         cancer = "",
         extra = unique(extra),
         extra2 = unique(extra2),
         dscore = 0)

vitamin.e <- vitamin.e %>%
  group_by(gene_symbol) %>%
  mutate(checked_gene_symbol = unique(checked_gene_symbol)) %>%
  ungroup() %>%
  group_by(checked_gene_symbol, show_drug_name) %>%
  mutate(pathways = unique(pathways),
         target_marker = unique(target_marker),
         ind_pathway = unique(ind_pathway),
         gene_dependency = unique(gene_dependency),
         gscore = unique(gscore),
         driver_gene = unique(driver_gene)) %>%
  ungroup() %>%
  group_by(checked_gene_symbol, show_drug_name, source, resistance) %>%
  mutate(alteration = unique(alteration)) %>%
  ungroup() %>%
  group_by(checked_gene_symbol, show_drug_name,source, alteration) %>%
  mutate(resistance = unique(resistance)) %>%
  ungroup() %>%
  unique()

pandrugs.new <- pandrugs.new %>%
  filter(show_drug_name != "VITAMIN E") %>%
  bind_rows(vitamin.e) %>%
  unique()

# Update MIDOSTAURIN, L-ASPARAGINASE and IDARUBICIN entries: cancer = blood
blood.drugs <- c("MIDOSTAURIN" = "TARGETED THERAPY", 
                 "L-ASPARAGINASE" = "CHEMOTHERAPY", 
                 "IDARUBICIN" = "CHEMOTHERAPY")

pandrugs.new <- pandrugs.new %>%
  mutate(is.blood = show_drug_name %in% names(blood.drugs),
         cancer = if_else(is.blood, "blood", cancer),
         extra2 = if_else(is.blood, blood.drugs[show_drug_name], extra2),
         dscore = case_when(is.blood & target_marker == "target" ~ 1,
                            is.blood & target_marker == "marker" ~ 0.9,
                            TRUE ~ dscore),
         dscore = if_else(is.blood & resistance == "resistance", -1 * dscore, 
                          dscore))

# Remove extra columns and check collapsed columns
pandrugs.new <- pandrugs.new %>%
  select(-c(version, is.blood)) %>%
  unique()

# Save table
write.table(pandrugs.new, col.names = TRUE, row.names = FALSE, quote = FALSE,
            sep = "\t", file = paste0(out.dir, "PanDrugs_Jun_06_2024_gs_dr.tsv"))

# Filter sensitivity/resistance entries
pandrugs.new <- pandrugs.new %>%
  filter(resistance != "sensitivity / resistance")

# Save table
write.table(pandrugs.new, col.names = TRUE, row.names = FALSE, quote = FALSE,
            sep = "\t", file = paste0(out.dir, "PanDrugs_Jun_06_2024_gs_dr_filtered.tsv"))

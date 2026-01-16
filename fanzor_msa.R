# Clear environment
rm(list = ls())
# Required packages
library(tidyverse)
library(ape)
# library(gridExtra)
# library(ggpubr)
if (require(showtext)) {
  showtext_auto()
  if (interactive()) showtext_opts(dpi = 100) else showtext_opts(dpi = 300)
}

# Set working directory
wd <- "/scratch1/kdeweese/jordan_fanzor/"
setwd(wd)

# Functions
# Checks existence of file or directory and errors if FALSE
checkPath <- function(test_path) {
  if (file.exists(test_path)) {
    return(TRUE)
  } else {
    stop(paste0("Cannot follow path (", test_path, ")."))
  }
}

# Input
# Species table file
species_tab_file <- "oyster_fanzor/species_table.tsv"
# Table of Fanzor orthologs (supplementary data)
# Bao, W. and Jurka, J. Mobile DNA 4, 12 (2013). https://doi.org/10.1186/1759-8753-4-12
ref_fnz_tab_file <- "Fz_TnpB_orthologs.tsv"
# Directory containing input protein FASTAs for OrthoFinder
prot_dir <- "proteomes"
# Path to results from OrthoFinder run
res_dirs <- list.files(
  path = paste(prot_dir, "OrthoFinder", sep = "/"),
  pattern = "Results_",
  full.names = T
)
# Use most recent OrthoFinder Results directory
res_dir <- file.info(res_dirs) %>%
  filter(isdir == "TRUE") %>%
  arrange(desc(mtime)) %>%
  slice_head(n = 1) %>%
  rownames()
# Main results table from OrthoFinder 
ortho_tsv_file <- paste(res_dir, "Orthogroups/Orthogroups.tsv", sep = "/")
checkPath(ortho_tsv_file)

# Import data
# Import species table
species_tab <- read_tsv(species_tab_file)
names(species_tab) <- gsub("# ", "", names(species_tab))
# Import reference Fanzor ortholog table
ref_fnz_tab <- read_tsv(ref_fnz_tab_file)
ref_fnz_tab <- ref_fnz_tab %>%
  separate(
    `start-end:strand`,
    into = c("start", "end", "strand"),
    sep = "[-:]",
    convert = T
  ) %>%
  mutate(strand = case_when(
    grepl("minus", strand) ~ "-",
    grepl("plus", strand) ~ "+"
  ))
# Import OrthoFinder results
ortho_tab <- read_tsv(ortho_tsv_file)
# Parse GFF file paths
species_tab <- species_tab %>%
  mutate(GFF = paste("ncbi_dataset/data", NCBI_Accession, "genomic.gff", sep = "/"))
fungus_gff_file <- species_tab %>%
  filter(Species == "Spizellomyces punctatus") %>%
  pull(GFF)
checkPath(fungus_gff_file)
fungus_gff <- read.gff(fungus_gff_file, GFF3 = T)
fungus_gff <- fungus_gff %>%
  mutate(across(where(is.factor), as.character)) %>%
  # Separate attributes col into columns named by regex (e.g., AC=## -> AC column)
  mutate(attributes = str_split(attributes, ";") %>%
           map(~ {
             key_vals <- str_split_fixed(.x, "=", 2)
             set_names(key_vals[, 2], key_vals[, 1])
           })) %>%
  unnest_wider(attributes, names_repair = "unique")

# Analysis
fungus_fz <- read_tsv("Fz_S_punctatus_prot_ids.txt", col_names = "protein_id")

fungus_gff_filt <- fungus_gff %>%
  # Filter for fungus orthologs
  filter(protein_id %in% fungus_fz)

ortho_tab %>%
  filter(grepl(paste(fungus_fz$protein_id, collapse = "|"), `Spizellomyces_punctatus-proteome`)) %>%
  View


# Clear environment
rm(list = ls())
# Required packages
library(ape)
library(Biostrings)
library(tidyverse)
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
    return(FALSE)
    # stop(paste0("Cannot follow path (", test_path, ")."))
  }
}
# Read protein FASTA file
readProtFa <- function(prot_fa) {
  prot <- readAAStringSet(prot_fa)
  prot_df <- tibble(
    name = names(prot),
    seq = as.character(prot),
    width = width(prot)
  )
  prot_df <- prot_df %>%
    mutate(
      Protein_ID = word(name, start = 1, end = 1),
      Product = gsub("^\\S* | \\[.*", "", name, perl = T),
      .before = name
    )
  return(prot_df)
}
# Extract protein IDs from GenBank file
getProts <- function(gb_file) {
  if (checkPath(gb_file)) {
    # Parse output filename
    out_fn <- gsub("\\.gb", "_prot_ids.tsv", gb_file)
    # Read output file if exists
    if (checkPath(out_fn)) {
      prot_ids <- read_tsv(out_fn)
    } else {
      gb <- read_lines(gb_file)
      prot_ids <- grep("protein_id=", gb, value = T)
      prot_ids <- gsub(".*protein_id=|\"", "", prot_ids)
      prot_ids <- tibble(Protein_ID = prot_ids)
      # Save parsed table to file
      write_tsv(x = prot_ids, file = out_fn)
    }
    return(prot_ids)
  }
}

# Input
# Species table file
species_tab_file <- "oyster_fanzor/species_table.tsv"
# # Table of Fanzor orthologs (supplementary data)
# # Bao, W. and Jurka, J. Mobile DNA 4, 12 (2013).
# # https://doi.org/10.1186/1759-8753-4-12
# ref_fnz_tab_file <- "Fz_TnpB_orthologs.tsv"
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
# # Import reference Fanzor ortholog table
# ref_fnz_tab <- read_tsv(ref_fnz_tab_file)
# ref_fnz_tab <- ref_fnz_tab %>%
#   separate(
#     `start-end:strand`,
#     into = c("start", "end", "strand"),
#     sep = "[-:]",
#     convert = T
#   ) %>%
#   mutate(strand = case_when(
#     grepl("minus", strand) ~ "-",
#     grepl("plus", strand) ~ "+"
#   ))
# Import species table
species_tab <- read_tsv(species_tab_file)
names(species_tab) <- gsub("# ", "", names(species_tab))
# Import OrthoFinder results
ortho_tab <- read_tsv(ortho_tsv_file)
# Import protein FASTAs
prot_df <- sapply(
  pull(species_tab, Proteome, Species),
  readProtFa,
  USE.NAMES = F,
  simplify = F
) %>%
  bind_rows(.id = "Species") # collapse into data frame
# Parse GFF file paths
species_tab <- species_tab %>%
  mutate(GFF = paste(
    "ncbi_dataset/data",
    NCBI_Accession,
    "genomic.gff",
    sep = "/"
  ))
# Parse Fanzor protein filenames
# Fanzor protein GenBank files from Saito, M., Nature, 620, 660–668 (2023).
# https://doi.org/10.1038/s41586-023-06356-2
species_tab <- species_tab %>%
  mutate(
    Fz_GenBank = paste0(
      "Fz_",
      paste(
        str_sub(Species, start = 1, end = 1),
        word(Species, start = 2, end = 2),
        sep = "_"
      ),
      ".gb"
    ),
    Fz_Proteins = paste0(
      "Fz_",
      paste(
        str_sub(Species, start = 1, end = 1),
        word(Species, start = 2, end = 2),
        sep = "_"
      ),
      "_prot_ids.tsv"
    )
  )
# Import Fanzor protein IDs from GenBank files
fz_prot_df <- sapply(
  pull(species_tab, Fz_GenBank, Species),
  getProts,
  USE.NAMES = F,
  simplify = F
) %>%
  bind_rows(.id = "Species") # collapse into data frame
# Keep only Fanzor proteins with matching IDs in proteomes
fz_prot_df <- fz_prot_df %>%
  filter(
    paste0(Species, Protein_ID) %in% paste0(prot_df$Species, prot_df$Protein_ID)
  )

# Analysis
# Separate out orthogroup protein IDs
ortho_long <- ortho_tab %>%
  pivot_longer(
    !Orthogroup,
    names_to = "Species",
    names_pattern = "(.*)-proteome",
    values_to = "Protein_ID"
  ) %>%
  separate_longer_delim(Protein_ID, delim = ", ") %>%
  mutate(Species = gsub("_", " ", Species))
# Annotate Fanzor orthologs
ortho_long <- ortho_long %>%
  mutate(Fanzor = case_when(
    paste0(Species, Protein_ID) %in%
      paste0(fz_prot_df$Species, fz_prot_df$Protein_ID) ~
      paste(Species, Protein_ID, sep = ";")
  )) %>%
  group_by(Orthogroup) %>%
  mutate(Fanzor_Ortho = paste(na.omit(Fanzor), collapse = ", ")) %>%
  ungroup()
# Filter for only Fanzor orthologs
ortho_fz <- ortho_long %>%
  filter(Fanzor_Ortho != "")

ortho_fz %>%
  distinct(Orthogroup)

ortho_fz %>%
  filter(if_all(Protein_ID, ~!is.na(.)), .by = c(Orthogroup, Species))

# fungus_gff_file <- species_tab %>%
#   filter(Species == "Spizellomyces punctatus") %>%
#   pull(GFF)
# checkPath(fungus_gff_file)
# fungus_gff <- read.gff(fungus_gff_file, GFF3 = T)
# fungus_gff <- fungus_gff %>%
#   mutate(across(where(is.factor), as.character)) %>%
#   # Separate attributes col into cols named by regex (e.g., AC=## -> AC column)
#   mutate(attributes = str_split(attributes, ";") %>%
#            map(~ {
#              key_vals <- str_split_fixed(.x, "=", 2)
#              set_names(key_vals[, 2], key_vals[, 1])
#            })) %>%
#   unnest_wider(attributes, names_repair = "unique")
# fungus_gff_filt <- fungus_gff %>%
#   # Filter for fungus orthologs
#   filter(protein_id %in% pull(
#     filter(fz_prot_df, Species == "Spizellomyces punctatus"),
#     protein_id
#   ))
# 



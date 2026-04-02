# Clear environment
rm(list = ls())
# Required packages
suppressPackageStartupMessages(
  library(tidyverse, quietly = T, warn.conflicts = F)
)
library(stringdist, quietly = T, warn.conflicts = F)
suppressPackageStartupMessages(library(msa, quietly = T, warn.conflicts = F))
library(Biostrings, quietly = T, warn.conflicts = F)
library(ape, quietly = T, warn.conflicts = F)
library(phangorn, quietly = T, warn.conflicts = F)
suppressPackageStartupMessages(library(ggtree, quietly = T, warn.conflicts = F))
suppressPackageStartupMessages(
  library(ggtreeExtra, quietly = T, warn.conflicts = F)
)
suppressPackageStartupMessages(library(ggmsa, quietly = T, warn.conflicts = F))
library(ggtext)
if (require(showtext, quietly = T, warn.conflicts = F)) {
  showtext_auto()
  if (interactive()) showtext_opts(dpi = 100) else showtext_opts(dpi = 300)
}

# Functions
# Checks existence of file or directory and errors if FALSE
checkPath <- function(test_path) {
  if (file.exists(test_path)) {
    return(TRUE)
  } else {
    stop(paste0("Cannot follow path (", test_path, ")."))
    return(FALSE)
  }
}
# Format species Latin name
formatSpc <- function(spc) {
  spc1 <- word(spc, 1, 1)
  spc1 <- substr(spc1, 1, 1)
  spc1 <- paste0(spc1, ".")
  spc2 <- word(spc, 2, 2)
  spc_f <- paste(spc1, spc2)
  return(spc_f)
}
# Get metadata from protein FASTA file
metaProtFa <- function(prot_fa) {
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
# Extract protein IDs from GenBank file (or existing output TSV file)
getProts <- function(gb_file) {
  if (file.exists(gb_file)) {
    gb <- read_lines(gb_file)
    prot_ids <- grep("protein_id=", gb, value = T)
    prot_ids <- gsub(".*protein_id=|\"", "", prot_ids)
    prot_ids <- tibble(Protein_ID = prot_ids)
    return(prot_ids)
  }
}
# Strip uninformative annotation parts
stripAnnots <- function(products) {
  products %>%
  str_replace_all("LOW QUALITY PROTEIN: ", "") %>%
    str_replace_all("\\-like", "") %>%
    str_replace_all(regex("isoform.*"), "") %>%
    str_replace_all("probable", "") %>%
    str_squish()
}
# Filter out completely uninformative annotations
filtAnnots <- function(products) {
  filt <- stripAnnots(products)
  filt[!str_detect(filt,
                       regex("^(uncharacterized|hypothetical)", ignore_case = T))]
}
# Find longest common prefix (word-boundary aware)
lcpLabel <- function(products) {
  inf <- filtAnnots(products)
  if (length(inf) == 0)
    return("uncharacterized protein")
  if (length(inf) == 1)
    return(inf)
  words <- str_split(inf, "\\s+")
  min_len <- min(lengths(words))
  common <- map_lgl(seq_len(min_len), function(i) {
    n_distinct(map_chr(words, i)) == 1
  })
  n_common <- which(!common)[1] - 1
  if (is.na(n_common))
    n_common <- min_len
  if (n_common == 0)
    return(inf[1])
  map_chr(words, ~ paste(.x[1:n_common], collapse = " ")) %>%
    unique() %>%
    dplyr::first()  # explicit namespace avoids conflicts
}
# Cluster products by string similarity, return cluster labels
clusterProds <- function(products, threshold = 0.2) {
  inf <- filtAnnots(products)
  if (length(inf) < 2)
    return(rep(1L, length(inf)))
  d <- stringdistmatrix(inf, method = "jw")
  # ward.D2 can produce non-monotone dendrograms on small/uniform matrices
  # fall back to "average" linkage which is more stable
  hc <- tryCatch(
    hclust(as.dist(d), method = "ward.D2"),
    error = function(e)
      hclust(as.dist(d), method = "average")
  )
  # cutree with 'h' fails if heights aren't strictly increasing (float precision issue)
  # sort the merge heights as a workaround
  hc$height <- sort(hc$height)
  cutree(hc, h = threshold)
}
# Returns one row per orthogroup with a combined label
collapseAnnots <- function(products, threshold = 0.2) {
  inf <- filtAnnots(products)
  n_total <- length(products)
  n_inf <- length(inf)
  # If nothing informative, just return that
  if (n_inf == 0)
    return("uncharacterized protein")
  # If only one unique informative name, use it
  if (n_distinct(inf) == 1)
    return(inf[1])
  clusters <- clusterProds(inf, threshold = threshold)
  # Label each cluster
  cluster_info <- tibble(product = inf, cluster = clusters) %>%
    summarize(n = n(), label = lcpLabel(product), .by = cluster) %>% 
    arrange(desc(n)) %>% # order clusters by size
    pull(label) %>%
    unique() %>% # create combined label
    paste(collapse = "; ")
}
# Run alignment
runMSA <- function(seqs) {
  # Strip terminal stop codons (*)
  seqs <- as(sub("\\*$", "", seqs), "AAStringSet")
  
  # Generate MSA object
  # Run MSA with protein sequences
  aln <- msa(seqs, method = "ClustalOmega")
  # Convert alignment to AAStringSet
  aln <- as(aln, "AAStringSet")
  
  # Generate tree object
  # Convert alignment to phyDat
  aln_phydat <- phyDat(as.matrix(aln), type = "AA")
  # Distance matrix
  dist <- dist.ml(aln_phydat)
  # Neighbor-Joining tree
  tree <- NJ(dist)
  
  return(list(aln, tree))
}
# Plot MSA alongside tree
plotTreeMSA <- function(prots, aln, tree, aln_title, label_offset, msa_width) {
  options(ignore.negative.edge = T) # ignore warnings about negative tree lengths
  # Extract metadata stored in prots (AAStringSet)
  meta <- data.frame(
    label = tree$tip.label, # must match tip labels in tree
    Species = mcols(prots)$Species,
    Product = mcols(prots)$Product
  ) %>%
    mutate(Product = gsub("'", "\\'", Product, fixed = T)) # escape quotes
  aln_sub <- paste0(
    "offset: ", label_offset, " (", max_chars, " char) ",
    "width: ", msa_width, " (", aln_length, " bp)"
  )
  tip_size <- 3.88  # ggtree default
  # Only reduce for very large trees where text would overlap vertically
  if (n_seqs > 40) tip_size <- max(2.5, 3.88 * (40 / n_seqs))
  # Adjust label offset by tip size
  label_offset <- (max_chars / char_per_off) * (tip_size / 3.88)
  p <- ggtree(tree) %<+% meta +
    geom_tiplab(aes(label = paste0(
      "italic('", formatSpc(Species), "')", " ~'", Product, "'"
    )), parse = T, align = T, size = tip_size)
  p_msatree <- msaplot(p, aln, offset = label_offset, width = msa_width) +
    labs(title = aln_title, subtitle = aln_sub) +
    theme(legend.position = "none")
  options(ignore.negative.edge = F) # turn warnings back on for user
  return(p_msatree)
}
# Round to nearest 10 (https://stackoverflow.com/a/25495249)
roundUp <- function(x) {
  if(length(x) != 1) stop("'x' must be of length 1")
  round(x + 5, digits = -1)
}
# Plot ggmsa sequence alignment with seqlogo and consensus bar graph
plotMSA <- function(prots, aln, aln_title) {
  # Add species names
  meta <- tibble(
    label = tree$tip.label, # must match tip labels in tree
    names = paste(mcols(prots)$Species, label)
    # names = paste(mcols(prots)$Species, mcols(prots)$Product)
  ) %>%
    pull(names, label)
  names(aln) <- meta[names(aln)]
  field <- roundUp(round(sqrt(max(width(aln)) * length(aln)) * 2, 0))
  # Generate ggmsa object
  p_msa <- ggmsa(aln, char_width = 0.7, seq_name = T, by_conservation = T) +
    labs(title = aln_title) +
    facet_msa(field = field)
    # geom_seqlogo() +
    # geom_msaBar()
  return(p_msa)
}

# Input
# Only take command line input if not running interactively
if (interactive()) {
  wd <- "/project2/noujdine_61/kdeweese/jordan_fanzor/"
  setwd(wd)
  # Species table file
  species_tab_file <- "oyster_fanzor/species_table.tsv"
  # URL list
  url_list_file <- "oyster_fanzor/url_list.txt"
  # Directory containing input protein FASTAs for OrthoFinder
  prot_dir <- "proteomes"
  # Output directory
  outdir <- "oyster_fanzor/"
} else {
  line_args <- commandArgs(trailingOnly = T)
  species_tab_file <- line_args[1]
  url_list_file <- line_args[2]
  prot_dir <- line_args[3]
  outdir <- line_args[4]
}
# Output
if (!dir.exists(outdir)) outdir <- "./" # if outdir does not exist, set to "./"
fz_ortho_tab_file <- "fanzor_ortho_table.tsv"
plot_dir <- "msa_plots"
fz_ortho_tab_file <- paste0(outdir, fz_ortho_tab_file)
plot_dir <- paste0(outdir, plot_dir)
dir.create(plot_dir, showWarnings = F) # create directory for plots, if needed

# Import datasets
# Import species table
checkPath(species_tab_file)
species_tab <- read_tsv(species_tab_file, show_col_types = F)
names(species_tab) <- names(species_tab) %>% # clean up header
  gsub("# ", "", .) %>%
  gsub(" ", "_", .)
# Import protein FASTAs
prot_df <- sapply(
  pull(species_tab, Proteome, Species),
  metaProtFa,
  USE.NAMES = F,
  simplify = F
) %>%
  bind_rows(.id = "Species") # collapse into data frame
# Import downloaded files
checkPath(url_list_file)
extra_files <- basename(read_lines(url_list_file))
extra_files <- paste0(outdir, extra_files) # find files in outdir
# Import table of clam gene annotation updates between reference versions
clam_spc <- "Mercenaria mercenaria"
clam_annot_corr_file <- grep(
  species_tab %>%
    filter(Species == clam_spc) %>%
    pull(Assembly_Accession),
  extra_files,
  value = T
)
checkPath(clam_annot_corr_file)
clam_annot_corr <- read_tsv(clam_annot_corr_file, show_col_types = F)
colnames(clam_annot_corr) <- gsub(" ", "_", colnames(clam_annot_corr))
# Remaining files in list are Fanzor protein GenBank files from supplement of
# Saito, M. et al. Nature, 620, 660–668, 2023. doi:10.1038/s41586-023-06356-2
fz_gbs <- grep(clam_annot_corr_file, extra_files, invert = T, value = T)
sapply(fz_gbs, checkPath)
gb_lines <- sapply(fz_gbs, read_lines, skip = 1, n_max = 1) # read 2nd line only
fz_gbs <- sapply(species_tab$Species, function(x) {
  y = names(grep(x, gb_lines, value = T))
  if (length(y) == 0) NA_character_ else y
})
species_tab <- species_tab %>%
  mutate(Fz_GenBank = unname(fz_gbs[Species])) # add column
# Import Fanzor protein IDs from GenBank files
fz_prot_list <- sapply(
  pull(species_tab, Fz_GenBank, Species),
  getProts,
  simplify = F
)
# Correct clam Fanzor protein IDs
clam_annot_corr <- clam_annot_corr %>%
  filter(
    previous_protein_accession %in% fz_prot_list[[clam_spc]]$Protein_ID
  ) %>%
  select(current_protein_accession, previous_protein_accession)
fz_prot_list[[clam_spc]] <- fz_prot_list[[clam_spc]] %>%
  left_join(
    clam_annot_corr, by = c("Protein_ID" = "previous_protein_accession")
  ) %>%
  select(current_protein_accession) %>%
  filter(!is.na(current_protein_accession))
colnames(fz_prot_list[[clam_spc]]) <- c("Protein_ID")
fz_prot_df <- fz_prot_list %>%
  # Collapse list into data frame
  bind_rows(.id = "Species") %>%
  # Keep only Fanzor proteins with matching IDs in proteome FASTAs
  filter(
    paste0(Species, Protein_ID) %in% paste0(prot_df$Species, prot_df$Protein_ID)
  )
# Import OrthoFinder results
res_dirs <- list.files(
  path = paste(prot_dir, "OrthoFinder", sep = "/"),
  pattern = "Results_",
  full.names = T
)
if (length(res_dirs) == 0) stop("No OrthoFinder results found.")
res_dir <- file.info(res_dirs) %>%
  filter(isdir == "TRUE") %>%
  # Use most recent Results dir
  arrange(desc(mtime)) %>%
  slice_head(n = 1) %>%
  rownames()
ortho_tsv_file <- paste(res_dir, "Orthogroups/Orthogroups.tsv", sep = "/")
checkPath(ortho_tsv_file)
ortho_tab <- read_tsv(ortho_tsv_file, show_col_types = F)

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
  # Mark specific proteins from Saito et al.
  mutate(Fanzor = case_when(
    paste0(Species, Protein_ID) %in%
      paste0(fz_prot_df$Species, fz_prot_df$Protein_ID) ~
      paste(Species, Protein_ID, sep = ";")
  )) %>%
  # Per orthogroup, Fanzor_Ortho = "Species1;Protein_ID1, Species2;Protein_ID2",
  # where each listed protein is Fanzor ID'd by Saito et al. found in Orthogroup
  group_by(Orthogroup) %>%
  mutate(Fanzor_Ortho = paste(na.omit(Fanzor), collapse = ", ")) %>%
  ungroup()
# Filter for only Fanzor orthologs in Magallana gigas
ortho_fz <- ortho_long %>%
  filter(Fanzor_Ortho != "") %>%
  group_by(Orthogroup) %>%
  filter(any(Species == "Magallana gigas" & !is.na(Protein_ID))) %>%
  ungroup() %>%
  mutate(
    Saito_Fanzor = case_when(!is.na(Fanzor) ~ "Yes", .default = "No"),
    .after = Fanzor
  ) %>%
  select(!Fanzor) %>%
  filter(!is.na(Protein_ID)) %>%
  # Add protein sequence and metadata
  left_join(prot_df, by = join_by(Species, Protein_ID)) %>%
  # Arrange for near-single-copy orthologs first
  mutate(n_Per_Orthogroup = n(), .by = Orthogroup) %>%
  mutate(n_Per_Species = n(), .by = c(Orthogroup, Species)) %>%
  arrange(n_Per_Orthogroup, desc(n_Per_Species)) %>%
  # Create consensus protein product labels for each orthogroup
  mutate(
    Ortho_Product = collapseAnnots(Product, threshold = 0.2),
    .by = Orthogroup,
    .after = Product
  )

# Write table of protein orthologs
write_tsv(x = ortho_fz, file = fz_ortho_tab_file)

# # Subset proteins for MSA
# prots_per_spc <- list()
# for (spc in species_tab$Species) {
#   prot_file <- species_tab %>%
#     filter(Species == spc) %>%
#     pull(Proteome)
#   filt_ids <- ortho_fz %>%
#     filter(Species == spc) %>%
#     pull(Protein_ID)
#   prots <- readAAStringSet(prot_file)
#   names(prots) <- word(names(prots), 1, 1)
#   prots <- prots[names(prots) %in% filt_ids]
#   names(prots) <- paste(formatSpc(spc), names(prots))
#   prots_per_spc[[spc]] <- prots
# }

# 
ortho_fz

# Group by orthogroup
ortho_split <- split(ortho_fz, ortho_fz$Orthogroup)
idxs <- which(names(ortho_split) %in% c("OG0000950", "OG0007788", "OG0000145", "OG0001120"))
# idx <- which(names(ortho_split) == "OG0007788")
# idx <- which(names(ortho_split) == "OG0000145")
# idx <- which(names(ortho_split) == "OG0001120")
idx <- which(names(ortho_split) == "OG0000950")
graph_stats <- list()
# for (og in names(ortho_split)[1:4]) {
# for (og in names(ortho_split)[idx]) {
# for (og in names(ortho_split)) {
for (og in names(ortho_split)[idxs]) {
  og_df <- ortho_split[[og]]
  og_prod <- og_df %>% distinct(Ortho_Product) %>% pull(Ortho_Product)
  prots <- og_df %>%
    select(seq, Protein_ID) %>%
    { AAStringSet(setNames(.$seq, .$Protein_ID)) }
  mcols(prots)$Species <- formatSpc(og_df$Species)
  mcols(prots)$Product <- og_df$Product
  res <- runMSA(prots)
  aln <- res[[1]]
  tree <- res[[2]]
  plot_title <- paste(og, og_prod, sep = ": ")
  
  # Set plot parameters based on input dataset
  # Maximum label character width
  max_chars <- max(nchar(paste(mcols(prots)$Species, mcols(prots)$Product)))
  aln_length <- unique(width(aln)) # MSA length in bp
  char_per_off <- 20 # characters per ggtree x-unit
  bp_per_width <- 650 # alignment bp per msaplot width unit
  # Convert to plot units
  # Offset: distance between tree and MSA (scale for alignments < 500bp)
  # label_offset <- max_chars / char_per_off
  label_offset <- (max_chars / char_per_off) * min(1, aln_length / 500)
  # msa_width <- aln_length / bp_per_width 
  # Width: ratio of MSA to tree width
  # msa_width <- max(aln_length / bp_per_width, 0.5) # min 0.5 units
  msa_width <- max(aln_length / bp_per_width, label_offset * 0.8) # offset ratio
  
  # Set ggsave width and height based on input dataset
  # Plot height driven purely by number of sequences
  n_seqs <- length(aln)
  h_mt <- n_seqs * 0.12 + 5 # 0.12in per seq + 5in margin
  # Plot width has additive components
  x_unit_to_in <- 0.3 # in per ggtree x-unit (tuned against known-good plot)
  tree_xmax <- max(ggtree(tree)$data$x, na.rm = T)
  # tree_in <- tree_xmax * x_unit_to_in * 7
  tree_in <- max(tree_xmax * x_unit_to_in * 7, 1 + n_seqs * 0.02)
  # tree_in <- max(tree_xmax * x_unit_to_in, 2) # min 2in
  # label_in <- label_offset * x_unit_to_in # calculate via offset
  label_in <- max_chars * 0.085 # directly from characters
  ref_msa_in <- 10.3 / 2.71  # inches per msaplot width unit
  msa_in <- max(msa_width * ref_msa_in, 3)  # min 3in
  # w_mt <- max(tree_in + label_in + msa_in, 8) # min 8in wide
  w_mt <- max(tree_in + label_in + msa_in, 10 + n_seqs * 0.15) # min set by n_seqs
  
  # Generate plots
  # p_msatree <- plotTreeMSA(prots, aln, tree, plot_title, label_offset, msa_width)
  p_msa <- plotMSA(prots, aln, plot_title)
  
  # Name plot files
  msa_tree_plot_file <- paste0(plot_dir, "/", og, "_MSA_tree.png")
  msa_plot_file <- paste0(plot_dir, "/", og, "_MSA.png")
  
  # Report parameters
  graph_stats[[og]] <- tibble(
    n_char = max_chars,
    n_seqs = n_seqs,
    aln_length = aln_length,
    tree_xmax = tree_xmax,
    label_offset = label_offset,
    msa_width = msa_width,
    tree_in = tree_in,
    label_in = label_in,
    msa_in = msa_in,
    w_mt = w_mt,
    h_mt = h_mt
  )
  
  # Save plots
  showtext_opts(dpi = 300)
  # ggsave(msa_tree_plot_file, p_msatree, dpi = 300, height = h_mt, width = w_mt)
  
  h_m <- 8 * max(1, (n_seqs / 22)) # increase height by number of sequences
  field <- roundUp(round(sqrt(max(width(aln)) * length(aln)) * 2, 0))
  w_m <- 16 * max(1, 0.1 * (field / 140)) # increase width by field size
  ggsave(msa_plot_file, p_msa, dpi = 300, height = h_m, width = w_m)
  showtext_opts(dpi = 100)
}
graph_stats <- graph_stats %>%
  bind_rows(.id = "Orthogroup")
graph_stats %>% print(n = 100)
# graph_stats %>% arrange(aln_length) %>% print(n=100)
graph_stats %>% arrange(tree_xmax) %>% print(n=100)
good <- c(
  "OG0000110",
  "OG0000311"
)
bad <-c(
  "OG0000950",
  "OG0006492"
)
graph_stats %>%
  filter(Orthogroup %in% c(good, bad)) %>%
  mutate(
    Graph = case_when(Orthogroup %in% good ~ "good", Orthogroup %in% bad ~ "bad")
  )

# fungus_gff_file <- species_tab %>%
#   filter(Species == "Spizellomyces punctatus") %>%
#   pull(GFF)
# checkPath(fungus_gff_file)
# fungus_gff <- read.gff(fungus_gff_file, GFF3 = T)
# fungus_gff <- fungus_gff %>%
#   mutate(across(where(is.factor), as.character)) %>%
#   # Separate attributes col into cols named by regex (e.g., AC=## -> AC col)
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



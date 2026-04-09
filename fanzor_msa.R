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
library(pwalign, quietly = T, warn.conflicts = F)
if (require(showtext, quietly = T, warn.conflicts = F)) {
  showtext_auto()
  if (interactive()) {
    my_dpi <- dev.size("px")[1]/dev.size("in")[1]
    showtext_opts(dpi = my_dpi)
  } else {
    showtext_opts(dpi = 300)
  }
}

# Functions
# Checks existence of file or directory and errors if FALSE
checkPath <- function(test_path) {
  sys_msg <- paste("Found file:", test_path)
  er_msg <- paste0("Cannot follow path (", test_path, ").")
  if (!file.exists(test_path)) stop(er_msg)
  return(sys_msg)
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
  if (length(products) == 0)
    return("uncharacterized protein")
  if (length(products) == 1)
    return(products)
  words <- str_split(products, "\\s+")
  min_len <- min(lengths(words))
  common <- map_lgl(seq_len(min_len), function(i) {
    n_distinct(map_chr(words, i)) == 1
  })
  n_common <- which(!common)[1] - 1
  if (is.na(n_common))
    n_common <- min_len
  if (n_common == 0)
    return(products[1])
  map_chr(words, ~ paste(.x[1:n_common], collapse = " ")) %>%
    unique() %>%
    dplyr::first()  # explicit namespace avoids conflicts
}
# Cluster products by string similarity, return cluster labels
clusterProds <- function(products, threshold = 0.2) {
  if (length(products) < 2)
    return(rep(1L, length(products)))
  d <- stringdistmatrix(products, method = "jw")
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
  
  # Generate tree object from MSA (Basic Protocol 5, doi.org/10.1002/cpbi.96)
  d <- stringDist(aln, method = "hamming")/width(aln)[1]
  tree <- bionj(d)
  
  return(list("msa" = aln, "tree" = tree))
}
# Plot MSA alongside tree
plotTreeMSA <- function(prots, aln, tree, plot_title, spc_names) {
  # Extract metadata stored in prots (AAStringSet)
  meta <- mcols(prots) %>%
    as.data.frame() %>%
    rownames_to_column("label") %>%
    arrange(factor(label, levels = tree$tip.label)) # match order of tree labels
  
  # Color species consistently
  sp_cols <- setNames(scales::hue_pal()(length(spc_names)), formatSpc(spc_names))
  
  # Resize title text, as needed
  title_chars <- nchar(plot_title) # title char. width
  title_size <- 13.2 # ggplot2 default
  if (title_chars > 80)  title_size <- 11 # ggplot2 default
  plot_title <- str_wrap(plot_title, 127, whitespace_only = F) # wrap long title
  
  # Resize tip label text, as needed
  n_seqs <- length(aln) # number of aligned sequences
  msa_height <- 0.8 # ggtree default
  lab_size <- 11  # ggplot2 default
  tip_size <- 2
  # Reduce MSA height for alignments with very few sequences
  if (n_seqs < 10) msa_height <- 0.5
  # Reduce font size for very large trees where text would overlap vertically
  if (n_seqs > 50) lab_size <- 7
  # Increase tip point size for small trees
  if (n_seqs < 10) tip_size <- 4
  
  # Tree plot
  p <- ggtree(tree) %<+% meta +
    geom_tiplab(
      aes(label = Product, color = Species),
      as_ylab = T, # display tip labels as y-axis label
      align = T # align nodes to right side with dotted line
    ) +
    geom_tippoint(size = tip_size) +
    theme(legend.text = element_text(face = "italic")) # italicize species names
  # Get x-axis range
  pb <- ggplot_build(p)
  tip_offset <- diff(pb$layout$panel_params[[1]]$x.range)*0.02
  # Rebuild tree plot with small tip offset
  p <- ggtree(tree) %<+% meta +
    geom_tiplab(
      aes(label = Product, color = Species),
      as_ylab = T, # display tip labels as y-axis label
      align = T # align nodes to right side with dotted line
    ) +
    geom_tippoint(
      aes(
        x = x + tip_offset,
        color = Species,
        shape = Saito_Fanzor
      ), size = tip_size
    ) +
    scale_color_manual(values = sp_cols) +
    scale_shape_manual(values = c(17, 8)) +
    theme(legend.text = element_text(face = "italic")) # italicize species names
  
  # Tree + MSA plot
  p_msatree <- msaplot(p, aln, height = msa_height, offset = 0.02) +
    labs(title = plot_title) +
    guides(
      fill = "none", # don't display MSA fill legend
      shape = "none" # don't display tip point shape legend
    ) +
    coord_cartesian(clip = "off") +
    theme(
      legend.location = "plot", # legend relative to entire plot area
      legend.position = "top",
      plot.title = element_text(size = title_size),
      axis.text.y.right = element_text(margin = margin(r = 10), size = lab_size),
      plot.margin = margin(l = 20)
    )
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
    labs(title = aln_title)
    # facet_msa(field = field)
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
  # Table of canonical Fz protein IDs from Saito et al. (2023)
  canon_fz_file <- "oyster_fanzor/canon_fz.tsv"
  # Directory containing input protein FASTAs for OrthoFinder
  prot_dir <- "proteomes"
  # Output directory
  outdir <- "oyster_fanzor/"
} else {
  line_args <- commandArgs(trailingOnly = T)
  species_tab_file <- line_args[1]
  url_list_file <- line_args[2]
  canon_fz_file <- line_args[3]
  prot_dir <- line_args[4]
  outdir <- line_args[5]
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
sapply(fz_gbs, checkPath, USE.NAMES = F)
gb_lines <- sapply(fz_gbs, read_lines, skip = 1, n_max = 1) # read 2nd line only
fz_gbs <- sapply(species_tab$Species, function(x) {
  y = names(grep(x, gb_lines, value = T))
  if (length(y) == 0) NA_character_ else y
})
species_tab <- species_tab %>%
  mutate(Fz_GenBank = unname(fz_gbs[Species])) # add column
# Import canonical Fanzor protein IDs per species from Saito et al. (2023)
checkPath(canon_fz_file)
canon_fz <- read_tsv(canon_fz_file, show_col_types = F)
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
  # Mark specific proteins from Saito et al. (2023)
  # Canonical Fanzor proteins
  left_join(., canon_fz) %>%
  # All Fanzor orthologs
  mutate(Fanzor = case_when(
    paste0(Species, Protein_ID) %in%
      paste0(fz_prot_df$Species, fz_prot_df$Protein_ID) ~
      paste(Species, Protein_ID, sep = ";")
  )) %>% 
  # Per orthogroup, Fanzor_Ortho = "Species1;Protein_ID1, Species2;Protein_ID2",
  # where each listed protein is Fanzor ID'd by Saito et al. found in Orthogroup
  group_by(Orthogroup) %>%
  mutate(Fanzor_Ortho = paste(na.omit(Fanzor), collapse = ", ")) %>%
  ungroup() %>%
  mutate(
    Saito_Fanzor = case_when(!is.na(Fanzor) ~ "Yes", .default = "No"),
    .after = Fanzor
  )
# Filter for only Fanzor orthologs in Magallana gigas
ortho_fz <- ortho_long %>%
  filter(Fanzor_Ortho != "") %>%
  group_by(Orthogroup) %>%
  filter(any(Species == "Magallana gigas" & !is.na(Protein_ID))) %>%
  ungroup() %>%
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
# Canonical Fz-containing orthogroups
ortho_canon <- ortho_long %>%
  group_by(Orthogroup) %>%
  filter(!all(is.na(Fz))) %>%
  ungroup() %>%
  select(!Fanzor) %>%
  filter(!is.na(Protein_ID)) %>%
  # Add protein sequence and metadata
  left_join(prot_df, by = join_by(Species, Protein_ID)) %>%
  # Create consensus protein product labels for each orthogroup
  mutate(
    Ortho_Product = collapseAnnots(Product, threshold = 0.2),
    .by = Orthogroup,
    .after = Product
  )

# Write table of protein orthologs
write_tsv(x = ortho_fz, file = fz_ortho_tab_file)

# Mark language clustering products by cluster #
ortho_fz <- ortho_fz %>%
  mutate(Cluster = clusterProds(stripAnnots(Product)), .by = Orthogroup) %>%
  mutate(
    Cluster_Product = lcpLabel(stripAnnots(Product)),
    .by = c(Orthogroup, Cluster),
    .after = Ortho_Product
  )

# Save species names
spc_names <- unique(ortho_fz$Species)
# Align and cluster sequences in each orthogroup
ortho_split <- split(ortho_fz, ortho_fz$Orthogroup) # split by orthogroup
canon_split <- split(ortho_canon, ortho_canon$Orthogroup) # split by orthogroup
# for (og in names(canon_split)) {
for (og in names(ortho_split)) {
  print(paste("Plotting alignments of:", og))
  og_df <- ortho_split[[og]] # subset orthogroup df
  # og_df <- canon_split[[og]] # subset orthogroup df
  prots <- og_df %>% # convert df to AAStringSet object
    pull(seq, Protein_ID) %>%
    AAStringSet()
  
  # Run protein MSA and generate tree from alignment
  res <- runMSA(prots)
  aln <- res[["msa"]]
  tree <- res[["tree"]]
  
  # # Label protein clusters from alignment
  # og_df <- aln %>%
  #   as("AAStringSet") %>%
  #   as.matrix() %>%
  #   phyDat(type = "AA") %>%
  #   dist.ml() %>%
  #   hclust(method = "complete") %>%
  #   cutree(h = 0.2) %>%
  #   tibble(Protein_ID = names(.), Tree_Cluster = .) %>%
  #   left_join(og_df, ., by = join_by(Protein_ID))
  
  # Add metadata to AAStringSet object
  mcols(prots)$Species <- formatSpc(og_df$Species) # formatted species name
  mcols(prots)$Product <- og_df$Product # protein product name
  mcols(prots)$Saito_Fanzor <- og_df$Saito_Fanzor # if protein is ID'd Fanzor
  # mcols(prots)$Cluster_Product <- og_df$Cluster_Product # cluster # of ortho name
  # mcols(prots)$Tree_Cluster <- og_df$Tree_Cluster # cluster # by alignment
  
  # Title plot with collapsed orthogroup product name
  og_prod <- og_df %>%
    distinct(Ortho_Product) %>%
    pull(Ortho_Product)
  plot_title <- paste(og, og_prod, sep = ": ") # plot title

  # Generate plots
  p_msatree <- plotTreeMSA(prots, aln, tree, plot_title, spc_names)
  # p_msa <- plotMSA(prots, aln, plot_title)

  # Name plot files
  msa_tree_plot_file <- paste0(plot_dir, "/", og, "_MSA_tree.png")
  msa_plot_file <- paste0(plot_dir, "/", og, "_MSA.png")
  
  # Save plots
  showtext_opts(dpi = 300)
  # Tree + MSA plot
  ggsave(msa_tree_plot_file, p_msatree, dpi = 300, height = 10, width = 10)
  # h_m <- 8 * max(1, (n_seqs / 22)) # increase height by number of sequences
  # field <- roundUp(round(sqrt(max(width(aln)) * length(aln)) * 2, 0))
  # w_m <- 16 * max(1, 0.1 * (field / 140)) # increase width by field size
  # ggsave(msa_plot_file, p_msa, dpi = 300, height = h_m, width = w_m)
  showtext_opts(dpi = my_dpi)
}

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



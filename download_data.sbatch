#!/bin/bash
#SBATCH -J download_data
#SBATCH -o %x_%j.log
#SBATCH --mem=10gb
#SBATCH --time=01:00:00

# Script variables
[[ -n "$SLURM_JOB_NAME" ]] || SLURM_JOB_NAME="download_data" # job name
scripts_dir="oyster_fanzor/" # scripts directory
species_list_file="species_list.txt" # input list of species
url_list_file="url_list.txt" # input list of URLs
prot_dir="proteomes" # output dir of proteome FASTAS
species_tab="species_table.tsv" # output species metadata table

# Redirect output if directory exists
[[ -d "$scripts_dir" ]] && species_tab="${scripts_dir}${species_tab}"

# Help message
if [[ "$1" = "-h" ]] || [[ "$1" = "--help" ]]
then
  echo "\
Takes as input a text file that lists one species Latin name per line.
Uses NCBI CLI's 'download datasets' to download current accessions for each
spceies, then unzips downloads and verifies file integrity. Copies proteome
FASTAs, renamed by species, to output directory $prot_dir for later OrthoFinder
analysis. Creates $species_tab, a TSV file with one species per line and
associated metadata with columns:
  # Species    Assembly Accession    Assembly Release Date    Proteome

Then, downloads necessary files from URLs listed in $url_list_file, including
GenBank files of Fanzor proteins from Saito, M. et al. (2023, doi:10.1038/s41586-023-06356-2)
and a comparison table for the current M. mercenaria genome annotations (GCF_021730395.1),
as the M. mercenaria GenBank file was built with the older genome (GCF_014805675.1).

Input:
  - species_list (default: $species_list_file)
  - url_list (default: $url_list_file)

Usage: sbatch ${SLURM_JOB_NAME}.sbatch <species_list> <url_list>

Output:
  ├── ncbi_dataset.zip
  ├── md5sum.txt
  └── ncbi_dataset/
      └── data/
          ├── assembly_data_report.jsonl
          ├── dataset_catalog.json
          └── <accession_IDs>/
          │   ├── genomic.gff
          │   ├── *_genomic.fna
          │   ├── protein.faa
          │   ├── rna.fna
          │   ├── cds_from_genomic.fna
          │   └── sequence_report.jsonl
  ├── $prot_dir/
      └── <species>_proteome.faa
  └── $scripts_dir (if exists, else: ./)
      ├── 41586_2023_6356_MOESM6_ESM.txt
      ├── 41586_2023_6356_MOESM8_ESM.txt
      ├── 41586_2023_6356_MOESM9_ESM.txt
      ├── 41586_2023_6356_MOESM10_ESM.txt
      ├── GCF_021730395.1_MADL_Memer_1_compare_prev.txt.gz
      └── $species_tab

Requires:
 - NCBI Datasets command-line tools (CLI) v18.14.0 (https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/)
 - UnZip v6.00+
 - curl v8+"
  exit 0
fi

# Print date and time
date
echo

# Input
# Species list file
if [[ -n "$1" ]]
then
  species_list_file="$1"
else
  [[ -d "$scripts_dir" ]] && species_list_file="${scripts_dir}${species_list_file}"
fi
# URL list file contents:
# 4 Fanzor protein GenBank files (Saito, M. et al. doi:10.1038/s41586-023-06356-2)
# 1 annotation comparison for M. mercenaria genome (GCF_021730395.1)
# Fanzor GenBank files built with older M. mercenaria genome (GCF_014805675.1)
if [[ -n "$2" ]]
then
  url_list_file="$2"
else
  [[ -d "$scripts_dir" ]] && url_list_file="${scripts_dir}${url_list_file}"
fi

# -e  exit on error
# -u  treat unset variables as errors
# -o pipefail  catch errors in pipes
set -euo pipefail

# Check input files exist
if [[ -f "$species_list_file" && -f "$url_list_file" ]]
then
  echo "Using species in list: $species_list_file"
  mapfile -t species_list < <(grep -v "#" "$species_list_file")
  echo "Found ${#species_list[@]} species:"
  printf '%s\n' "${species_list[@]}"
  echo
  echo "Using URLs in: $url_list_file"
  mapfile -t url_list < "$url_list_file"
  echo "Found ${#url_list[@]} URLs:"
  printf '%s\n' "${url_list[@]}"
else
  echo "Error: at least one expected input ($species_list_file or $url_list_file) not found."
  exit 1
fi
echo

# Load conda
cond=~/.conda_for_sbatch.sh
if [[ -f "$cond" ]]
then
  source "$cond"
else
  echo "Error on source of $cond"
  exit 1
fi

# Load NCBI CLI
module purge
conda activate ncbi
datasets --version || exit 1

# Search NCBI for reference accession numbers for each species
# Specify fields to query from NCBI
fields_arr=("accession" "assminfo-release-date")
fields="$(IFS=,; echo "${fields_arr[*]}")"
# Create species table from output
echo "Creating species table ($species_tab)..."
for i in "${!species_list[@]}"
do
  spc="${species_list[i]}"
  if [[ "$i" -eq 0 ]]
  then
    datasets summary genome taxon "$spc" --annotated --as-json-lines | dataformat tsv genome > test.tsv
    datasets summary genome taxon "$spc" --reference --annotated --as-json-lines \
      | dataformat tsv genome --fields "$fields" \
      > "$species_tab"
  else
    datasets summary genome taxon "$spc" --reference --annotated --as-json-lines \
      | dataformat tsv genome --fields "$fields" --elide-header \
      >> "$species_tab"
  fi
done
print_spc=("# Species" "${species_list[@]}") # add header to array
paste <(printf '%s\n' "${print_spc[@]}") "$species_tab" > tmp.tsv
mv tmp.tsv "$species_tab"
# Array of accession numbers
mapfile -t access_list < <(grep -v "#" "$species_tab" | cut -f2)
for i in "${!species_list[@]}"; do grep "${access_list[i]}" "${scripts_dir}"*ESM.txt; done
for i in "${!species_list[@]}"; do echo "${access_list[i]}"; done
# Download accessions
echo "Downloading accessions from NCBI..."
set -x # echo commands
datasets download genome accession "${access_list[@]}" \
  --include gff3,rna,cds,protein,genome,seq-report \
  --no-progressbar
set +x
conda deactivate

# Load zip module
module load unzip
unzip -v || exit 1

# Unzip NCBI dataset
set -x # echo commands
unzip -o ncbi_dataset.zip -x README.md
set +x
echo

# Verify integrity after unzip
if md5sum -c md5sum.txt
then
  echo
  echo "Download and decompression of NCBI accessions successful."
else
  echo
  echo "Error: Download and decompression of NCBI accessions unsuccessful."
  exit 1
fi
echo

# Create directory to copy proteomes to
mkdir -p "$prot_dir"
# printf "# %s\t%s\t%s\n" "Species" "NCBI_Accession" "Proteome" > "${species_tab}_tmp"
prot_list=() # array of paths to be column 'Proteome' in species_tab
# Copy each proteome to directory and rename by species
for i in "${!species_list[@]}"
do
  species="${species_list[i]}"
  access="${access_list[i]}"
  species_fn="${species// /_}" # remove spaces from species name
  og_proteome="ncbi_dataset/data/$access/protein.faa"
  proteome="$prot_dir/${species_fn}-proteome.faa"
  if [[ -f "$og_proteome" ]]
  then
    set -x
    cp "$og_proteome" "$proteome"
    set +x
    # Save new proteome path to array
    prot_list+=("$proteome")
  else
    echo "Error: proteome FASTA ($og_proteome) not found for $species ($access)."
  fi
done
echo
# Add header
print_prots=("Proteome" "${prot_list[@]}")
# Check for expected number of paths in array
if [[ "${#prot_list[@]}" -eq "${#species_list[@]}" ]]
then
  # Update table with new column
  new_col="Proteome"
  echo "Updating table ($species_tab) with $new_col column..."
  # If column already exists, remove and replace
  if grep -q "^#.*$new_col" "$species_tab"
  then
    # Find existing column index
    col_idx="$(grep "#" "$species_tab" \
      | tr '\t' '\n' \
      | grep -n "$new_col" \
      | cut -d: -f1)"
    # Remove old column and add new one
    paste \
      <(cut -f "$col_idx" --complement "$species_tab") \
      <(printf '%s\n' "${print_prots[@]}") \
      > tmp.tsv
  else
    # Append new column to TSV with paste
    paste "$species_tab" <(printf '%s\n' "${print_prots[@]}") > tmp.tsv
  fi
  mv tmp.tsv "$species_tab"
else
  echo "Error: check FASTAs in $prot_dir/"
  # exit 1
fi
echo

# Check curl
module purge
module load curl
curl --version || exit 1

# Download files from URLs with curl
for i in "${!url_list[@]}"
do
  url="${url_list[i]}"
  if [[ -d "$scripts_dir" ]]
  then
    curl -O --output-dir "$scripts_dir" "$url"
  else
    curl -O "$url"
  fi
done

# Print date and time
echo "Finished at: $(date)"

#!/bin/bash
#SBATCH -J download_data
#SBATCH --mem=10gb
#SBATCH --time=01:00:00
#SBATCH -o %x.log

# Script variables
prot_dir="proteomes"

# Help message
if [[ "$1" = "-h" ]] || [[ "$1" = "--help" ]]
then
  [[ -n "$SLURM_JOB_NAME" ]] || SLURM_JOB_NAME="download_data"
  echo "\
Takes as input two-column TSV file (species_tab), where column 1 contains 
species Latin names and column 2 has corresponding NCBI accession IDs.
Script downloads with NCBI CLI 'download datasets', unzips download, and 
verifies integrity of files. Then, copies proteome FASTAs, renamed by species,
to output directory $prot_dir, and updates <species_tab> to include new paths.

Usage: sbatch ${SLURM_JOB_NAME}.sbatch <species_tab>

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
  └── proteomes/
      └── <species>_proteome.faa

Requires:
 - NCBI Datasets command-line tools (CLI) v18.14.0 (https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/)
 - zip v3.0+"
  exit 0
fi

# Print date and time
date_fmt="%-I:%M:%S %p (%a %d %b %Y)" # date format
date +"$date_fmt"
echo

# Input
# Species table
species_tab="oyster_fanzor/species_table.tsv"
if [[ -f "$species_tab" ]]
then
  echo "Using species table file: $species_tab"
  mapfile -t species_list < <(grep -v "#" "$species_tab" | awk -F '\t' '{print $1}')
  mapfile -t access_list < <(grep -v "#" "$species_tab" | awk -F '\t' '{print $2}')
  echo "Found ${#species_list[@]} species: ${species_list[*]}"
  echo "Found ${#access_list[@]} accessions: ${access_list[*]}"
else
  echo "Error: species table file ($species_tab) not found."
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

# Load NCBI CLI conda env
module purge
conda activate ncbi
datasets --version

# Download accessions
cmd=(
  datasets
  download
  genome
  accession
  "${access_list[@]}"
  --include gff3,rna,cds,protein,genome,seq-report
  --no-progressbar
)
echo "${cmd[*]}"
"${cmd[@]}"

# Unzip NCBI dataset
cmd2=(unzip -o ncbi_dataset.zip -x README.md)
echo "${cmd2[*]}"
"${cmd2[@]}"

# Verify integrity after unzip
if md5sum -c md5sum.txt
then
  echo
  echo "Download and decompression of NCBI accessions successful."
  echo
else
  echo
  echo "Error - download and decompression of NCBI accessions unsuccessful."
  exit 1
fi

# Copy proteomes to directory and rename by species
mkdir -p "$prot_dir"
# Wipe previous tempfile and add column header
printf "# %s\t%s\t%s\n" "Species" "NCBI_Accession" "Proteome" > "${species_tab}_tmp"
for i in "${!species_list[@]}"
do
  species="${species_list[i]}"
  access="${access_list[i]}"
  species_fn="${species// /_}" # remove spaces from species name
  og_proteome="ncbi_dataset/data/$access/protein.faa"
  proteome="$prot_dir/${species_fn}-proteome.faa"
  if [[ -f "$og_proteome" ]]
  then
    cmd3=(cp "$og_proteome" "$proteome")
    echo "${cmd3[*]}"
    "${cmd3[@]}"
    printf "%s\t%s\t%s\n" "$species" "$access" "$proteome" >> "${species_tab}_tmp"
  fi
done

# Update species table with new relative paths
og_lines="$(wc -l < "$species_tab")"
new_lines="$(wc -l < "${species_tab}_tmp")"
if [[ "$og_lines" -eq "$new_lines" ]]
then
  echo
  echo "Updating species table with proteome paths:"
  cmd4=(mv "${species_tab}_tmp" "$species_tab")
  echo "${cmd4[*]}"
  "${cmd4[@]}"
  echo
else
  echo
  echo "Error - check ${species_tab}_tmp."
  echo
fi

# Print date and time
date_fmt="%-I:%M:%S %p (%a %d %b %Y)" # date format
date +"$date_fmt"
echo

#!/bin/bash
#SBATCH -J fanzor_msa
#SBATCH -o %x_%j.log
#SBATCH --mem=50gb
#SBATCH --time=05:00:00
#SBATCH --cpus-per-task=12

# Script variables
[[ -n "$SLURM_JOB_NAME" ]] || SLURM_JOB_NAME="fanzor_msa" # job name
r_script="${SLURM_JOB_NAME}.R" # R script (named from job)
scripts_dir="oyster_fanzor/" # scripts directory
species_tab="species_table.tsv" # input species metadata table (default filename)
url_list_file="url_list.txt" # input list of URLs
canon_fz_file="canon_fz.tsv" # input table of canonical Fz protein IDs
prot_dir="proteomes" # input dir of proteome FASTAS
if [[ -d "$scripts_dir" ]] # redirect if scripts dir exists
then
  r_script="${scripts_dir}${r_script}"
  species_tab="${scripts_dir}${species_tab}"
  url_list_file="${scripts_dir}${url_list_file}"
  canon_fz_file="${scripts_dir}${canon_fz_file}"
fi

# Help message
if [[ "$1" = "-h" ]] || [[ "$1" = "--help" ]]
then
  echo "\
Takes as input table of species metadata and aligns Fanzor orthologs found using
OrthoFinder.

Usage: sbatch ${SLURM_JOB_NAME}.sbatch <species_tab>

Requires:
 - curl (https://curl.se/)
 - R v4.5.2 (https://www.r-project.org/)"
  exit 0
fi

# Print date and time
date
echo

# Input
# Species table
[[ -n "$1" ]] && species_tab="$1"

# -e  exit on error
# -u  treat unset variables as errors
# -o pipefail  catch errors in pipes
set -euo pipefail

# Check input file exists
if [[ -f "$species_tab" ]]
then
  echo "Using species table file: $species_tab"
  # Array of proteome paths
  mapfile -t prot_list < <(grep -v "#" "$species_tab" | cut -f4) 
  echo "Found ${#prot_list[@]} protein FASTAs."
else
  echo "Error: species table file ($species_tab) not found."
  exit 1
fi

# Load R modules
module purge
module load apptainer/1.4.5 rstats/4.5.2
Rscript --version || exit 1

# Run MSA of Fanzor orthologs in R
set -x
Rscript "$r_script" "$species_tab" "$url_list_file" "$canon_fz_file" \
  "$prot_dir" "$scripts_dir"
set +x

# Print date and time
echo "Finished at: $(date)"

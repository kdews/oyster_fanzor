#!/bin/bash
#SBATCH -J orthofinder
#SBATCH -o %x_%j.log
#SBATCH --mem=25gb
#SBATCH --time=05:00:00
#SBATCH --cpus-per-task=32

# Script variables
[[ -n "$SLURM_JOB_NAME" ]] || SLURM_JOB_NAME="orthofinder" # job nmw

# Help message
if [[ "$1" = "-h" ]] || [[ "$1" = "--help" ]]
then
  echo "\
Takes as input directory of protein FASTAs intended for OrthoFinder.

Usage: sbatch ${SLURM_JOB_NAME}.sbatch <species_tab>

Requires:
 - OrthoFinder v3.1.2 (https://github.com/OrthoFinder/OrthoFinder)"
  exit 0
fi

# Print date and time
date_fmt="%-I:%M:%S %p (%a %d %b %Y)" # date format
date +"$date_fmt"
echo

# Parallelization
# Low-mem threads
thr_t="$SLURM_CPUS_PER_TASK"
# High-mem threads
thr_a="$(( SLURM_CPUS_PER_TASK / 8 ))"

# Input
# Proteome diretctory
if [[ -n "$1" ]]
then
  prot_dir="$1"
else
  prot_dir="proteomes"
fi

# -e  exit on error
# -u  treat unset variables as errors
# -o pipefail  catch errors in pipes
set -euo pipefail

# Check directory exists
if [[ -d "$prot_dir" ]]
then
  echo "Using proteome diretctory: $prot_dir"
else
  echo "Error: proteome diretctory ($prot_dir) not found."
  exit 1
fi

# Load conda
module purge
cond=~/.conda_for_sbatch.sh
if [[ -f "$cond" ]]
then
  source "$cond"
else
  echo "Error on source of $cond"
  exit 1
fi

# Check existing directory for log of completed OrthoFinder run
echo
outdir="${prot_dir}/OrthoFinder"
[[ -d "$outdir" ]] && mapfile -t log < <(ls -t "${outdir}"/Results*/Log.txt)
for i in "${!log[@]}"
do
  if [[ -f "${log[i]}" ]] && grep -q "OrthoFinder run completed" "${log[i]}"
  then
    echo "Found previous successful run: ${log[i]}"
    res_dir="${log[i]%%"Log.txt"}"
    echo "Skipping OrthoFinder step, results in: $res_dir"
    break
  fi
done
# Run OrthoFinder
if [[ -z "$res_dir" ]]
then
  # Load OrthoFinder conda env
  echo
  conda activate orthofinder3
  orthofinder -v
  cmd="orthofinder -t $thr_t -a $thr_a -f $prot_dir"
  echo "$cmd"
  $cmd
  mapfile -t res_list < <(ls -t "${outdir}"/Results*/)
  res_dir="${res_list[0]}"
  conda deactivate
  echo "OrthoFinder results: $res_dir"
fi

# Print date and time
echo "Finished at: $(date)"

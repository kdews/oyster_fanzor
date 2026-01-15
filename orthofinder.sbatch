#!/bin/bash
#SBATCH -J orthofinder
#SBATCH --mem=25gb
#SBATCH --time=05:00:00
#SBATCH --cpus-per-task=32
#SBATCH -o %x.log

# Script variables
prot_dir="proteomes"
# R_mod_file="modules_for_Rscript.txt"
# r_script="${SLURM_JOB_NAME}.R"

# Help message
if [[ "$1" = "-h" ]] || [[ "$1" = "--help" ]]
then
  [[ -n "$SLURM_JOB_NAME" ]] || SLURM_JOB_NAME="orthofinder"
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

# # Prepend scripts_dir
# R_mod_file="${scripts_dir}${R_mod_file}"
# r_script="${scripts_dir}${r_script}"

# Parallelization
# Low-mem threads
thr_t="$SLURM_CPUS_PER_TASK"
# High-mem threads
thr_a="$(( SLURM_CPUS_PER_TASK / 8 ))"

# Input
# Species table
species_tab="oyster_fanzor/species_table.tsv"
if [[ -f "$species_tab" ]]
then
  echo "Using species table file: $species_tab"
  mapfile -t prot_list < <(grep -v "#" "$species_tab" | awk -F '\t' '{print $3}')
  echo "Found ${#prot_list[@]} protein FASTAs."
else
  echo "Error: species table file ($species_tab) not found."
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
    echo "Skipping OrthoFinder step, using results in: $res_dir"
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

# # Load R modules
# module purge
# IFS=" " read -r -a R_mods <<< "$(cat "$R_mod_file")"
# module load "${R_mods[@]}"

# # Run analysis of OrthoFinder results in R
# cmd="Rscript $r_script $species_tab $prot_dir/ $res_dir Saccharina_latissima $scripts_dir"
# echo "$cmd"
# # $cmd

# Print date and time
date_fmt="%-I:%M:%S %p (%a %d %b %Y)" # date format
date +"$date_fmt"
echo

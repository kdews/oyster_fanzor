#!/bin/bash

mapfile -t files < <(ls *.gb)
for file in "${files[@]}"
do
  file_base="${file%%.gb}"
  new_file="${file_base}_prot_ids.txt"
  grep "protein_id=" "$file" | sed 's/.*protein_id=//g; s/"//g' > "$new_file"
done

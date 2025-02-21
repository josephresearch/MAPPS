#!/bin/bash

if [[ -z "$1" ]]; then
    echo "Usage: $0 <drug_name>"
    exit 1
fi

drug=$1

mkdir -p params

# Read sigma values from logP.txt into a temporary file
sigma_file="sigma_values.tmp"
awk '{if (NR > 1) print $1, $3}' ../Step1-Get_SM-SM_Parameters/logP.txt > "$sigma_file"

file="../Step2-Get-SM-Y-Parameters/2.2-Fit_CG_to_AA_6Y/data/CG_best_fit/sm_ff_${drug}.dat"

if [[ -f "$file" ]]; then
    eps_0=$(awk '/pair_coeff  9  41/ {print $5}' "$file")

    if [[ -n "$eps_0" ]]; then
        sigma=$(grep -i "^$drug " "$sigma_file" | awk '{print $2}')
        if [[ -n "$sigma" ]]; then
            drug_ff="params/${drug}_ff.txt"
            python get_contacts_final.py --drug "$drug" --eps_0 "$eps_0" --sigma "$sigma" > "$drug_ff"
        else
            echo "Warning: No sigma value found for $drug in logP.txt"
        fi
    else
        echo "Warning: No eps_0 value found for $drug in $file"
    fi
else
    echo "Warning: File not found for $drug: $file"
fi

# Clean up temporary file
rm -f "$sigma_file"

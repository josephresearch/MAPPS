#!/usr/bin/env bash

base_dir="./"
gmx_dir="$base_dir/gmx"

echo "Current working directory: $base_dir"

mkdir -p "$gmx_dir"

for chain_ind in {2..8}; do
    mkdir -p "$gmx_dir/seq$chain_ind"
    for i in {1..3}; do
        dest_dir="$gmx_dir/seq${chain_ind}/seq${chain_ind}_$i"
        src_dir="$base_dir/seq${chain_ind}/seq${chain_ind}_$i"

        echo "Creating directories in: $dest_dir"
        mkdir -p "$dest_dir/num_contacts"
        mkdir -p "$dest_dir/min_dist"
        
        echo "Checking source directory: $src_dir"
        if [ -d "$src_dir" ]; then
            echo "Source directory $src_dir exists. Moving files..."
            mv "$src_dir"/numcont* "$dest_dir/num_contacts" 2>/dev/null
            if [ $? -ne 0 ]; then
                echo "No numcont* files found in $src_dir."
            fi
            mv "$src_dir"/*.xvg "$dest_dir/min_dist" 2>/dev/null
            if [ $? -ne 0 ]; then
                echo "No *.xvg files found in $src_dir."
            fi
        else
            echo "Directory $src_dir does not exist. Skipping."
        fi
    done
done

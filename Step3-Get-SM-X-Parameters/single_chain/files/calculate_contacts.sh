#!/usr/bin/env bash

folders=("./")
 
for folder in "${folders[@]}"; do
    for chain_ind in {2..8}; do

        for i in {1..3}; do

            cp get_contacts_per_res.sh $folder/seq$chain_ind/seq${chain_ind}_$i
            cp make_res_ndx.sh $folder/seq$chain_ind/seq${chain_ind}_$i
            cp mkndx_0.sh $folder/seq$chain_ind 
            cp mkndx.sh $folder/seq$chain_ind 
        done

        cd "$folder/seq$chain_ind/"
        pwd

        bash mkndx_0.sh
        bash mkndx.sh
        
        cd ../
    done
done

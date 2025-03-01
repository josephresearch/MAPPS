#!/usr/bin/env bash

folders=('2-methylpropane')

for folder in "${folders[@]}"; do
    mkdir -p "$folder"
    cp files/* "$folder"
    
    for file in "$folder"/*; do
        sed -i "s/drug/$folder/g" "$file"
    done

    for chain_ind in {2..8}; do
        for i in {1..3}; do
            mkdir -p "$folder/seq$chain_ind/seq${chain_ind}_$i"
            cp "seq_$chain_ind/packmol.gro" "$folder/seq$chain_ind/seq${chain_ind}_$i"
            cp "seq_$chain_ind/chain_20_$chain_ind.top" "$folder/seq$chain_ind/seq${chain_ind}_$i/chain_20_$chain_ind.top"
            cp "seq_$chain_ind/chain_20_$chain_ind.top" "$folder/seq$chain_ind/seq${chain_ind}_$i/chain_20_0.top"
            cp "seq_$chain_ind/chain_20_$chain_ind.gro" "$folder/seq$chain_ind/seq${chain_ind}_$i"

            cp "$folder/topol.top" "$folder/seq$chain_ind/seq${chain_ind}_$i"
            cp "$folder/run_sim.sh" "$folder/seq$chain_ind/seq${chain_ind}_$i"
            cp "$folder/solvate.sh" "$folder/seq$chain_ind/seq${chain_ind}_$i"
            cp "$folder/npt.mdp" "$folder/seq$chain_ind/seq${chain_ind}_$i"
            cp "$folder/nvt.mdp" "$folder/seq$chain_ind/seq${chain_ind}_$i"
            cp "$folder/em.mdp" "$folder/seq$chain_ind/seq${chain_ind}_$i"
            
            cp "2-methylpropane/seq$chain_ind/seq${chain_ind}_1/npt.gro" "$folder/seq$chain_ind/seq${chain_ind}_$i"
            
            cd  "$folder/seq$chain_ind/seq${chain_ind}_$i"

            bash solvate.sh
            sbatch run_sim.sh
            cd ../../../
        done
    done
done

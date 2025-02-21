#!/usr/bin/env bash

for d in ./*/ ; do
    cd $d 
    sbatch run_nvt.sh
    cd ../
done

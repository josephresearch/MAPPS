#!/usr/bin/env bash

for d in ./seq*/ ; do
    cd $d
    bash make_res_ndx.sh
    cd ../
done


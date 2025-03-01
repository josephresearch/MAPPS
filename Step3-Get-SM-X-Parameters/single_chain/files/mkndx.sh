#!/usr/bin/env bash

for d in ./seq*/ ; do
    cd $d
    bash get_contacts_per_res.sh
    cd ../
done


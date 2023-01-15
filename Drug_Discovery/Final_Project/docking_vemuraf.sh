#!/usr/bin/env bash
#File:Docking


## Docking in 3og7 (chain A)
smina -r ./files/3og7_clean_prot.pdb -l ./files/docking_vemuraf_ligands.mol2 --center_x "1.868515" --center_y "-2.6376667" --center_z "-19.917727" --size_x "20" --size_y "20" --size_z "20" --num_modes 3 --exhaustiveness 8 --seed 8763585 -o ./files/docking_results_vemuraf_3og7.sdf > ./files/docking_rmsd_vemuraf_3og7.txt


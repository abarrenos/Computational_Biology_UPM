#!/usr/bin/env bash
#File:Docking

## Docking in 4xv2 (chain A and B)
smina -r ./files/4xv2_clean_prot.pdb -l ./files/docking_ligands.mol2 --center_x "-0.3668" --center_y "-3.1906571" --center_z "-20.677143" --size_x "20" --size_y "20" --size_z "20" --num_modes 3 --exhaustiveness 8 --seed 8763585 -o ./files/docking_results_4xv2_chA.sdf > ./files/docking_rmsd_4xv2_ChA.txt
smina -r ./files/4xv2_clean_prot.pdb -l ./files/docking_ligands.mol2 --center_x "-3.8276572" --center_y "-1.0898571" --center_z "6.133" --size_x "20" --size_y "20" --size_z "20" --num_modes 3 --exhaustiveness 8 --seed 8763585 -o ./files/docking_results_4xv2_chB.sdf > ./files/docking_rmsd_4xv2_ChB.txt 
## Docking in 3og7 (chain A)
smina -r ./files/3og7_clean_prot.pdb -l ./files/docking_ligands.mol2 --center_x "1.868515" --center_y "-2.6376667" --center_z "-19.917727" --size_x "20" --size_y "20" --size_z "20" --num_modes 3 --exhaustiveness 8 --seed 8763585 -o ./files/docking_results_3og7.sdf > ./files/docking_rmsd_3og7.txt


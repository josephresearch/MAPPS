gmx editconf -f packmol.pdb -o packmol.gro

gmx editconf -f packmol.gro -o  packmol.gro -c -box 4.0 4.0 4.0
gmx solvate -cp packmol.gro -cs ../../../../../../amber03ws.ff/tip4p2005.gro -o packmol_w.gro -p topol.top
echo q | gmx make_ndx -f packmol_w.gro

gmx grompp -f minim.mdp -c packmol_w.gro -r packmol_w.gro -p topol.top -o ions.tpr -maxwarn 10
echo SOL | gmx genion -s ions.tpr -p topol.top -conc 0.150 -neutral -o solvated.gro
gmx_gpu insert-molecules -f npt.gro -ci ../../../../small_molecules/drug/drug.gro -o molecule.gro -nmol 1 # -ip positions.dat

gmx_gpu solvate -cp molecule.gro -cs ../../../amber03ws.ff/tip4p2005.gro -o packmol_w.gro -p topol.top
echo q | gmx_gpu make_ndx -f packmol_w.gro

gmx_gpu grompp -f em.mdp -c packmol_w.gro -r packmol_w.gro -p topol.top -o ions.tpr -maxwarn 3
echo SOL | gmx_gpu genion -s ions.tpr -p topol.top -conc 0.150 -neutral -o solvated.gro
#!/bin/bash
#SBATCH --job-name=aa-chain           # create a short name for your job
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=8        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem=1G                # memory per node (4G per cpu-core is default)
#SBATCH --time=12:00:00          # total run time limit (HH:MM:SS)
#SBATCH --gres=gpu:1             # number of gpus per node

module purge
module load cudatoolkit/11.4
module load openmpi/gcc/4.1.0

gmx_gpu grompp -f em.mdp -c solvated.gro -r solvated.gro -p topol.top -o em.tpr -maxwarn 4
srun gmx_gpu mdrun -ntomp $SLURM_CPUS_PER_TASK -s em.tpr -v -deffnm em

gmx_gpu grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr -maxwarn 4
srun gmx_gpu mdrun -ntomp $SLURM_CPUS_PER_TASK -s nvt.tpr -v -deffnm nvt

gmx_gpu grompp -f npt.mdp -c em.gro -r em.gro -p topol.top -o npt.tpr -maxwarn 2
srun gmx_gpu mdrun -ntomp $SLURM_CPUS_PER_TASK -s npt.tpr -v -deffnm npt

gmx_gpu convert-tpr -s nvt.tpr -extend 1000000 -o nvt2.tpr
srun gmx_gpu mdrun -ntomp $SLURM_CPUS_PER_TASK -deffnm nvt2 -cpi nvt.cpt -noappend

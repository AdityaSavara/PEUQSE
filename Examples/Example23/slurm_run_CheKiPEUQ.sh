#!/bin/sh
##SBATCH --partition=general-compute --qos=general-compute
##SBATCH --partition=debug
##SBATCH --partition=largemem
#SBATCH --partition=debug --qos=debug
##SBATCH --partition=skylake --qos=skylake
#SBATCH --cluster=ub-hpc
##SBATCH --cluster=chemistry
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
##SBATCH --constraint=IB
#SBATCH --mem=128000
##SBATCH --mem=23000
#SBATCH --job-name=ex12neg1std_scaling2E3
#SBATCH --output=job.out
#SBATCH --mail-user=ericwalk@buffalo.edu
#SBATCH --mail-type=ALL
##SBATCH --requeue
module load python/py37-anaconda-2019.10
############
module list

python runfile_Example23_pos_1std_offset.py

#!/bin/bash
#SBATCH -J "MIX_MOD"
#SBATCH --output=/data/vkarthik/ortho/logs/slurm.%j.out
#SBATCH -N 1
#SBATCH -n 8

 ~/anaconda3/lib/R/bin/Rscript mixture_model.rscript

#!/bin/bash
#SBATCH --job-name=1D_model
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --time=24:00:00

# cargar entorno
module purge
module load gcc

# cargar entorno
source ~/.bashrc
conda activate guachos

# ejecutar
./1D_test.exe
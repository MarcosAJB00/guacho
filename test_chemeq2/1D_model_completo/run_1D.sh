#!/bin/bash
#SBATCH --job-name=1D_model
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --time=24:00:00

# cargar entorno
source ~/.bashrc
conda activate guacho

# ejecutar
./1D_test.exe
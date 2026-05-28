#!/bin/bash
#SBATCH --job-name=multimodel
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --time=24:00:00

# cargar entorno
module purge
module load gcc

# cargar entorno
source ~/.bashrc
conda activate guacho13

# ejecutar
./1D_test_multimodel.exe
#!/bin/bash

### Las líneas #SBATCH configuran los recursos de la tar
### (aunque parezcan estar comentadas)
### Nombre de la tarea

#SBATCH --job-name=hetau0
### Cola de trabajos a la cual enviar.
#SBATCH --ntasks=1
#SBATCH --partition=batch
###SBATCH --nodes=1
###SBATCH --ntasks-per-node=1
###SBATCH --cpus-per-task=1
###SBATCH --nodelist=clemente05,clemente06

### Tiempo de ejecucion. Formato dias-horas:minutos.
#SBATCH --time 7-00:00:00
#####24:00:00

### Script que se ejecuta al arrancar el trabajo
### Cargar el entorno del usuario incluyendo la funcionalidad de modules
. /etc/profile
### Cargar los módulos para la tare

###source /home/mbaracchi/.miniconda3/etc/profile.d/conda.sh

###conda activate guacho

module load gcc

module load openmpi

export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
###export OMPI_MCA_pml=ob1
###export OMPI_MCA_btl=self,vader,tcp
###export OMPI_MCA_btl=self,openib

srun he_lines_tau0
#srun he_lines_tau1
#srun he_lines_tau2

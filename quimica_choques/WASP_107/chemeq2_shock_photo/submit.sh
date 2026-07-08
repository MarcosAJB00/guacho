#!/bin/bash
#SBATCH --job-name=1D_photo        # nombre del job en la cola
#SBATCH --output=logs/run_%j.out      # stdout (%j = job ID)
#SBATCH --error=logs/run_%j.err       # stderr
#SBATCH --nodes=1                     # UN solo nodo (obligatorio para OpenMP)
#SBATCH --ntasks=1                    # un solo proceso MPI (no usamos MPI)
#SBATCH --cpus-per-task=1            # 5 cores para OpenMP
#SBATCH --time=12:00:00               # tiempo maximo (HH:MM:SS), ajustar
#SBATCH --partition=batch           # nombre de la particion, depende del cluster

# =======================================================================
# Configuracion de OpenMP
# =======================================================================
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK   # toma el valor de --cpus-per-task
export OMP_PROC_BIND=close                     # cada thread en un core cercano
export OMP_PLACES=cores                        # un thread por core fisico

# =======================================================================
# Informacion del job (util para debuggear)
# =======================================================================
echo "========================================"
echo "Job ID:       $SLURM_JOB_ID"
echo "Job name:     $SLURM_JOB_NAME"
echo "Nodo:         $SLURMD_NODENAME"
echo "Cores:        $OMP_NUM_THREADS"
echo "Inicio:       $(date)"
echo "========================================"

# =======================================================================
# Crear directorios de salida si no existen
# =======================================================================
mkdir -p output
mkdir -p logs

# =======================================================================
# Correr el ejecutable
# =======================================================================
./1D_test_multi_omp.exe

echo "========================================"
echo "Fin: $(date)"
echo "========================================"
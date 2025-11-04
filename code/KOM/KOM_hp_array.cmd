#!/bin/bash
#SBATCH -J kom_hp_array
#SBATCH --partition=ncpu
#SBATCH --array=1-50
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH --mem=5G

# send Slurm’s own batch stdout/err to nowhere (we’ll manage logs ourselves)
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

# === set the shared output directory ONCE ===
# Use an absolute, unique path per submission (job ID)
OUTDIR="${PWD}/kom_hp_${SLURM_ARRAY_JOB_ID}"
mkdir -p "${OUTDIR}"
mkdir -p "${OUTDIR}/out/"
mkdir -p "${OUTDIR}/err/"
mkdir -p "${OUTDIR}/res/"
# make OUTDIR visible to the R process
export OUTDIR
# math threads per task
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
# total number of iteration
export N_REPEAT=50


# run the R script inside the container, redirecting logs into OUTDIR
srun --cpu-bind=cores --exclusive \
  singularity exec ../../R-container.sif \
  Rscript KOM_hp_selec_script.R --outdir "$OUTDIR/res" \
  1>  "${OUTDIR}/out/${SLURM_ARRAY_TASK_ID}.log" \
  2>  "${OUTDIR}/err/${SLURM_ARRAY_TASK_ID}.log"

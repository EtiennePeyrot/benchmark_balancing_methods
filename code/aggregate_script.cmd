#!/bin/bash
#SBATCH -J aggregate_res
#SBATCH --partition=ncpu
#SBATCH --ntasks=1
#SBATCH --mem=50G

SIM="${1}"

# set the directory
SIMDIR="${PWD}/${SIM}"
RESDIR="${SIMDIR}/res"
export RESDIR


# run the R script inside the container
srun --cpu-bind=cores --exclusive \
  singularity exec ../R-container.sif \
  Rscript aggregate_script.R


echo "aggregation done!"
echo "begin compression to zip"

zip -9 -m -j "${SIM}.zip" "data.RData"

echo "end compression to zip"


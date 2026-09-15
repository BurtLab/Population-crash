#!/bin/bash
#PBS -l walltime=04:00:00
#PBS -l select=1:ncpus=2:mem=4gb
#PBS -j oe
#PBS -N subsample_tajimasD

module load miniforge/3
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate msprime_env

# -----------------------------------------------------------------------
# EDIT THIS: point at the directory containing the already-simulated
# .trees files for the seasonal model, 50% suppression, sigma=0.5 scenario.
# Adjust SCENARIO_GLOB to match your actual directory naming for the 50%
# crash level.
# -----------------------------------------------------------------------
INPUT_ROOT="$HOME/project/hetero"
SCENARIO_GLOB="$INPUT_ROOT/Seasonal_0.5/DryTime_*/hetero_0.5"

OUTPUT_DIR="$INPUT_ROOT/stat_hetero/Seasonal_0.5_subsample"
mkdir -p "$OUTPUT_DIR"

for SCENARIO_DIR in $SCENARIO_GLOB; do

  cd "$SCENARIO_DIR"

  # 1. Compute Tajima's D at n=10 and n=30, subsampled from the existing
  #    n=50 tree sequences, for every seed/model .trees file in this directory.
  python "$HOME/subsample_tajimasD.py" 10,30 1

  cd -

  # 2. Wrangle the per-seed subsample outputs into one combined CSV
  #    (Time, N, Seed, Model, Value) in OUTPUT_DIR.
  Rscript --vanilla -e "
    library(reticulate);
    use_python('$(which python)', required = TRUE);
    source('$HOME/wrangle_subsample.R');
    process_subsample_tajimasD('$SCENARIO_DIR', '$OUTPUT_DIR')
  "

done

# Move log files to logs/ directory
mv "$PBS_O_WORKDIR"/*.o* "$HOME/logs/"

#!/bin/bash
#PBS -l walltime=12:00:00
#PBS -l select=1:ncpus=2:mem=2gb
#PBS -j oe
#PBS -N wrangle

module load anaconda3/personal

# Activate msprime environment
source activate msprime_env

# Set working directories
mkdir -p $HOME/project/hetero/stat_hetero
INPUT_DIR="$HOME/project/hetero"
OUTPUT_DIR="$INPUT_DIR/stat_hetero"

# 1. Process Constant directories
for Ne_DIR in "$INPUT_DIR"/Constant_*; do
    for SIGMA_DIR in "$Ne_DIR"/hetero_*; do
        
        # Output path: stat/Constant_*/hetero_*
        NE_NAME=$(basename "$Ne_DIR")
        SIGMA_NAME=$(basename "$SIGMA_DIR")
        TARGET_DIR="$OUTPUT_DIR/$NE_NAME/$SIGMA_NAME"
        
        Rscript --vanilla -e "
          library(reticulate);
          use_python('$(which python)', required = TRUE);
          source('$HOME/wrangle.R'); 
          process_directory_stats('$SIGMA_DIR', '$TARGET_DIR')
        "
    done
done
   
# 2. Process Seasonal directories (Nested DryTime directories)
for Seasonal_DIR in "$INPUT_DIR"/Seasonal_*; do
  for DryTime_DIR in "$Seasonal_DIR"/DryTime_*; do
    for SIGMA_DIR in "$DryTime_DIR"/hetero_*; do

      # Output path: stat/Seasonal_*/DryTime_*/hetero_*
      SEASON_NAME=$(basename "$Seasonal_DIR")
      DRY_NAME=$(basename "$DryTime_DIR")
      SIGMA_NAME=$(basename "$SIGMA_DIR")
      TARGET_DIR="$OUTPUT_DIR/$SEASON_NAME/$DRY_NAME/$SIGMA_NAME"

      Rscript --vanilla -e "
        library(reticulate);
        use_python('$(which python)', required = TRUE);
        source('$HOME/wrangle.R');
        process_directory_stats('$SIGMA_DIR', '$TARGET_DIR')
      "
    done
  done
done
    

# Move log files to logs/ directory
mv "$PBS_O_WORKDIR"/*.o* "$HOME/logs/"
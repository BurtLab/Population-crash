#!/bin/bash
#PBS -l walltime=12:00:00
#PBS -l select=1:ncpus=2:mem=10gb
#PBS -J 1-200
#PBS -j oe
#PBS -N simulation

module load anaconda3/personal

# Activate msprime environment
source activate msprime_env

# Define Ne values after crash
half_period=(6)
crash_value=(0.1)
sigma_value=(0 0.1 0.3 0.5)

# Copy scripts to $TMPDIR (HPC temporary directory)
cp $HOME/{seasonal.py,decline_seasonal.py,statistics.py} $TMPDIR

# Run the simulations for each Ne value
for half_period in "${half_period[@]}"; do
    for crash in "${crash_value[@]}"; do
        for sigma in "${sigma_value[@]}"; do
            python seasonal_hetero.py $half_period $sigma
            python decline_seasonal_hetero.py $half_period $sigma $crash
            python statistics.py

            mkdir -p "$HOME/project/hetero/Seasonal_$crash/DryTime_$half_period/hetero_$sigma"
            mv "$TMPDIR/control_"* "$TMPDIR/intervention_"* "$HOME/project/hetero/Seasonal_$crash/DryTime_$half_period/hetero_$sigma/"
        done
    done
done


# Move log files to logs directory
mkdir -p $HOME/logs
mv "$HOME/"*.o* $HOME/logs/
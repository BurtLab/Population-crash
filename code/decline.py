import msprime
import math
import os
import sys
import random

seed_number = int(os.getenv("PBS_ARRAY_INDEX"))  # Find out the job number


### Parameters ###
M = 30000
sigma = float(sys.argv[1])
logM = math.log(M) - 0.5 * sigma**2
crash_percent = float(sys.argv[2])
Ne = int(round(random.lognormvariate(logM, sigma)))
Ne_current = int(Ne * crash_percent)  # Population size after crash
seqlength = 1e6
mu = 2.5e-8
recombination_rate = 10 * mu 
timepoints = list(range(0, 49, 3))

### Define demographic model ###
demography = msprime.Demography()
demography.add_population(initial_size=Ne_current)  # The most recent population size
demography.add_population_parameters_change(
    time=36,  # Time of population decline (36 generations ago)
    initial_size=Ne,  # Population size before crash
    population=0
)

samples = [msprime.SampleSet(50, population=0,time=t) for t in timepoints]

### Define map positions and rates ###
r_chrom = recombination_rate
r_break = math.log(2)
chrom_positions = [0, seqlength, seqlength*2]
map_positions = [
    chrom_positions[0],
    chrom_positions[1],
    chrom_positions[1] + 1,
    chrom_positions[2]
]
rates = [r_chrom, r_break, r_chrom]
rate_map = msprime.RateMap(position=map_positions, rate=rates)


### Simulate ancestry ###
ts = msprime.sim_ancestry(
    samples=samples,
    demography=demography,  # Use the defined demographic model
    sequence_length=seqlength*2,
    recombination_rate=rate_map,
    model=[
        msprime.DiscreteTimeWrightFisher(duration=50),
        msprime.StandardCoalescent()
        ],
    random_seed=seed_number
)

# Add mutations
mts = msprime.sim_mutations(
    ts,
    rate=mu,
    random_seed=seed_number
)


### Save tree sequence ###
mts.dump(f"intervention_{seed_number}.trees")


### Record Ne ###
with open(f"intervention_Ne_{seed_number}.txt", "a") as f:
    f.write(f"{Ne}")
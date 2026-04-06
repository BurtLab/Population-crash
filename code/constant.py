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
Ne = int(round(random.lognormvariate(logM, sigma)))
seqlength = 1e6  # Per chromosome
mu = 2.5e-8  # Mutation rate (theta=4Ne*mu, theta=0.03)
recombination_rate = 10 * mu
timepoints = list(range(0, 49, 3))  # Time points from 0 to 48 in steps of 3
samples = [msprime.SampleSet(50, time=t) for t in timepoints]

### Define map positions and rates ###
r_chrom = recombination_rate
r_break = math.log(2)  # High recombination at chromosome boundary (~independent 0.5)
chrom_positions = [0, seqlength, seqlength*2] # Chromosome endpoints
map_positions = [
    chrom_positions[0],     # 0
    chrom_positions[1],     # 1e6
    chrom_positions[1] + 1, # 1e6+1 (gap)
    chrom_positions[2]      # 2e6
]
rates = [r_chrom, r_break, r_chrom]
rate_map = msprime.RateMap(position=map_positions, rate=rates)


### Simulate ancestry ###
ts = msprime.sim_ancestry(
    samples=samples,
    population_size=Ne,
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
mts.dump(f"control_{seed_number}.trees")


### Record Ne ###
with open(f"control_Ne_{seed_number}.txt", "a") as f:
    f.write(f"{Ne}")
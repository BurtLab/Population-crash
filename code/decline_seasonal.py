import msprime
import math
import os
import sys
import random

half_period = int(sys.argv[1])
seed_number = int(os.getenv("PBS_ARRAY_INDEX"))


### Parameters ###
M = 30000
sigma = float(sys.argv[2])
logM = math.log(M) - 0.5 * sigma**2
crash_percent = float(sys.argv[3])
Ne_wet = int(round(random.lognormvariate(logM, sigma)))
Ne_dry = int(Ne_wet * 0.1)
max_time = 1e6 
seqlength = 1e6
mu = 2.5e-8
recombination_rate = 10 * mu
timepoints = list(range(0, 49, 3))


### Define demographic model ###
demography = msprime.Demography()
demography.add_population(initial_size=Ne_wet*0.1)
sim_time = 0
sim_Ne = Ne_wet

while sim_time <= max_time:
    # Scale Ne by crash_percent for times before 36
    scale = crash_percent if sim_time < 36 else 1.0
    
    demography.add_population_parameters_change(
        time=sim_time,
        initial_size=sim_Ne * scale,
        population=0,
    )
    
    # Switch Ne
    sim_Ne = Ne_wet if sim_Ne == Ne_dry else Ne_dry
    sim_time += half_period

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
    ts,  # Input tree sequence
    rate=mu,  # Mutation rate
    random_seed=seed_number
)


### Save tree sequence ###
mts.dump(f"intervention_{seed_number}.trees")

### Record Ne ###
with open(f"intervention_Ne_{seed_number}.txt", "a") as f:
    f.write(f"{Ne_wet}")
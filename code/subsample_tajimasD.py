"""
subsample_tajimasD.py

Computes Tajima's D at reduced sample sizes (e.g. n=10, 30 individuals per
cluster) by randomly subsampling individuals from the existing n=50 tree
sequences already simulated for a given scenario - no new msprime runs
needed (sampling consistency of the coalescent: a random subset of
individuals from a larger simulated sample has the same distribution as an
independent simulation at that smaller sample size).

Run from within a directory containing the *.trees files for one scenario.

Usage:
    python subsample_tajimasD.py 10,30 1
    (first arg: comma-separated n_values, second arg: RNG seed)
"""

import sys
import glob

import numpy as np
import pandas as pd
import tskit

# Time points from 0 to 48 in steps of 3 - matches statistics.py's `timepoints`
TIMEPOINTS = list(range(0, 49, 3))


def get_tree_prefixes():
    return [f.replace(".trees", "") for f in glob.glob("*.trees")]


def subsample_individual_nodes(ts, t, n_sub, rng):
    sample_ids = ts.samples(time=t)
    individual_ids = np.unique([ts.node(nid).individual for nid in sample_ids])
    chosen = rng.choice(individual_ids, size=n_sub, replace=False)
    sub_nodes = np.concatenate([ts.individual(i).nodes for i in chosen])
    return sub_nodes


def tajimas_d_for_subsample(ts, t, n_sub, rng):
    nodes = subsample_individual_nodes(ts, t, n_sub, rng)
    sub_ts = ts.simplify(samples=nodes)
    return sub_ts.Tajimas_D()


def process_file(prefix, n_values, rng):
    ts = tskit.load(f"{prefix}.trees")

    records = []
    for t in TIMEPOINTS:
        for n_sub in n_values:
            d = tajimas_d_for_subsample(ts, t, n_sub, rng)
            records.append({"Time": t, "N": n_sub, "Tajimas_D": d})

    df = pd.DataFrame(records)
    df.to_csv(f"{prefix}_subsample_tajimasD.csv", index=False)


def main():
    n_values = [int(x) for x in sys.argv[1].split(",")]
    seed = int(sys.argv[2])

    rng = np.random.default_rng(seed)

    for prefix in get_tree_prefixes():
        process_file(prefix, n_values, rng)


if __name__ == "__main__":
    main()

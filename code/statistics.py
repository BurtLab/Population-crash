import tskit
import pandas as pd
import numpy as np
import glob


# Time points from 0 to 72 in steps of 3
timepoints = list(range(0, 49, 3))

# Get the prefix of the tree file
def get_tree_prefixes():
    prefixes = []
    for f in glob.glob("*.trees"):
        prefixes.append(f.replace(".trees", ""))
    return prefixes

prefix = get_tree_prefixes()


##### Calculate most summary statistics #####

def oneway_stats(ts):
    stats = {
        "Number of trees": ts.num_trees,
        "Number of mutations": ts.num_mutations,
        "Density of segregating sites": ts.segregating_sites(),
        "Nucleotide diversity (pi)": ts.diversity(),
        "Tajima's D": ts.Tajimas_D(),
    }
    return stats

def oneway_output(ts, prefix):
    results = {}
    for t in timepoints:
        # Extract samples for the current time point
        sample_ids = ts.samples(time=t)
        sub_ts = ts.simplify(samples=sample_ids)

        stats = oneway_stats(sub_ts)
        results[t] = stats

    oneway_df = pd.DataFrame(results).T  # Transpose to have time points as rows
    oneway_df.index.name = "t"
    oneway_df.to_csv(f"{prefix}.csv")



##### LD calculation functions (C-Optimized) #####
import ctypes
import os
import subprocess

# ---------------------------------------------------------
# 1. C Source Code
# ---------------------------------------------------------
C_CODE = """
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct {
    double mean;
    double median;
    double std;
    double min;
    double max;
    long long n_pairs;
} Stats;

int compare(const void * a, const void * b) {
    float fa = *(const float*)a;
    float fb = *(const float*)b;
    return (fa > fb) - (fa < fb);
}

void compute_ld_stats(
    const char *geno1, int n_sites1,
    const char *geno2, int n_sites2,
    int n_ind,
    Stats *out_stats
) {
    long long n_pairs_total = (long long)n_sites1 * n_sites2;
    if (n_pairs_total == 0) return;
    
    // Allocate site stats arrays
    double *dp1 = malloc(n_sites1 * sizeof(double));
    double *dh1 = malloc(n_sites1 * sizeof(double));
    double *denom1 = malloc(n_sites1 * sizeof(double));
    
    double *dp2 = malloc(n_sites2 * sizeof(double));
    double *dh2 = malloc(n_sites2 * sizeof(double));
    double *denom2 = malloc(n_sites2 * sizeof(double));
    
    if(!dp1 || !dh1 || !denom1 || !dp2 || !dh2 || !denom2) return;

    // Precompute Set 1
    for(int i=0; i<n_sites1; ++i) {
        int sum_p = 0;
        int sum_h = 0;
        for(int k=0; k<n_ind; ++k) {
            char val = geno1[i*n_ind + k];
            sum_p += val;
            if(val == 2) sum_h++;
        }
        dp1[i] = (double)sum_p / (2.0 * n_ind);
        dh1[i] = (double)sum_h / n_ind;
        denom1[i] = dp1[i] - 2.0*dp1[i]*dp1[i] + dh1[i];
    }
    
    // Precompute Set 2
    for(int j=0; j<n_sites2; ++j) {
        int sum_p = 0;
        int sum_h = 0;
        for(int k=0; k<n_ind; ++k) {
            char val = geno2[j*n_ind + k];
            sum_p += val;
            if(val == 2) sum_h++;
        }
        dp2[j] = (double)sum_p / (2.0 * n_ind);
        dh2[j] = (double)sum_h / n_ind;
        denom2[j] = dp2[j] - 2.0*dp2[j]*dp2[j] + dh2[j];
    }
    
    // Allocate results for median calculation
    float *results = malloc(n_pairs_total * sizeof(float));
    long long count = 0;
    
    double sum = 0.0;
    double sum_sq = 0.0;
    float min_val = 0.0f; 
    float max_val = 0.0f;
    int first = 1;

    // Main Loop
    for(int i=0; i<n_sites1; ++i) {
        if(denom1[i] <= 1e-12) continue; 
        
        for(int j=0; j<n_sites2; ++j) {
            if(denom2[j] <= 1e-12) continue;

            int delta = 0;
            const char *g1_row = &geno1[i*n_ind];
            const char *g2_row = &geno2[j*n_ind];
            
            for(int k=0; k<n_ind; ++k) {
                delta += g1_row[k] * g2_row[k];
            }
            
            double ddelta = (double)delta / (2.0 * n_ind);
            ddelta = ddelta - 2.0 * dp1[i] * dp2[j];
            ddelta = ddelta * (double)n_ind / (n_ind - 1.0);
            
            double denom = denom1[i] * denom2[j];
            float r2 = (float)((ddelta * ddelta) / denom);
            
            results[count++] = r2;
            
            sum += r2;
            sum_sq += r2*r2;
            if (first) { min_val=r2; max_val=r2; first=0; }
            else {
                if(r2 < min_val) min_val = r2;
                if(r2 > max_val) max_val = r2;
            }
        }
    }
    
    out_stats->n_pairs = count;
    if(count > 0) {
        out_stats->mean = sum / count;
        out_stats->std = sqrt(sum_sq/count - (sum/count)*(sum/count));
        out_stats->min = min_val;
        out_stats->max = max_val;
        
        qsort(results, count, sizeof(float), compare);
        if (count % 2 == 0) {
            out_stats->median = (results[count/2 - 1] + results[count/2]) / 2.0;
        } else {
            out_stats->median = results[count/2];
        }
    } else {
        out_stats->mean = 0; out_stats->std = 0;
        out_stats->min = 0; out_stats->max = 0;
        out_stats->median = 0;
    }

    free(dp1); free(dh1); free(denom1);
    free(dp2); free(dh2); free(denom2);
    free(results);
}
"""

# ---------------------------------------------------------
# 2. Python Parts
# ---------------------------------------------------------

def get_unphased_genotypes(ts, sample_indices):
    """
    Extract unphased genotypes from tree sequence for specified samples.
    """
    geno_matrix = ts.genotype_matrix()
    haplotypes = geno_matrix[:, sample_indices]
    n_individuals = len(sample_indices) // 2
    n_sites = haplotypes.shape[0]
    
    # 0, 1, 2 genotype
    unphased = np.zeros((n_sites, n_individuals), dtype=np.int8)
    for i in range(n_individuals):
        unphased[:, i] = haplotypes[:, 2*i] + haplotypes[:, 2*i + 1]
    
    return unphased

def calculate_maf_vectorized(genotypes):
    """
    Vectorized MAF calculation for filtering.
    """
    n_samples = genotypes.shape[1]
    freq = np.sum(genotypes, axis=1) / (2 * n_samples)
    maf = np.minimum(freq, 1 - freq)
    return maf

# ---------------------------------------------------------
# 3. Compilation and Wrapper
# ---------------------------------------------------------

class Stats(ctypes.Structure):
    _fields_ = [("mean", ctypes.c_double),
                ("median", ctypes.c_double),
                ("std", ctypes.c_double),
                ("min", ctypes.c_double),
                ("max", ctypes.c_double),
                ("n_pairs", ctypes.c_longlong)]

def compile_and_load_clib():
    lib_name = "./libld_fast.so"
    if not os.path.exists("ld_fast.c"):
        with open("ld_fast.c", "w") as f:
            f.write(C_CODE)
            
    if not os.path.exists(lib_name):
        cmd = ["gcc", "-O3", "-shared", "-fPIC", "-o", lib_name, "ld_fast.c"]
        subprocess.check_call(cmd)

    lib = ctypes.CDLL(lib_name)
    lib.compute_ld_stats.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.int8, flags='C_CONTIGUOUS'), ctypes.c_int,
        np.ctypeslib.ndpointer(dtype=np.int8, flags='C_CONTIGUOUS'), ctypes.c_int,
        ctypes.c_int,
        ctypes.POINTER(Stats)
    ]
    return lib

# Initialize Library
try:
    _LD_LIB = compile_and_load_clib()
except Exception as e:
    print(f"Warning: Could not compile C library: {e}")
    _LD_LIB = None

# ---------------------------------------------------------
# 4. Main Function
# ---------------------------------------------------------

def calculate_ld_stats(ts, sample_id_dict, maf_threshold=0.05):

    seqlength = 1e6
    all_sites = np.arange(ts.num_sites)
    chrom1_mask = ts.sites_position < seqlength
    chrom2_mask = ts.sites_position > seqlength
    
    chrom1_sites = all_sites[chrom1_mask]
    chrom2_sites = all_sites[chrom2_mask]
    
    summary_results = []
    
    # Iterate over timepoints (keys of dictionary)
    for t in sample_id_dict.keys():
        sample_indices = sample_id_dict[t]
        
        # 1. Get Genotypes
        unphased_geno = get_unphased_genotypes(ts, sample_indices)
        n_individuals = unphased_geno.shape[1]
        
        # 2. Filter Sites
        maf_all = calculate_maf_vectorized(unphased_geno)
        c1_indices = chrom1_sites[maf_all[chrom1_sites] > maf_threshold]
        c2_indices = chrom2_sites[maf_all[chrom2_sites] > maf_threshold]
        
        # 3. Prepare C-Contiguous Arrays
        geno1 = np.ascontiguousarray(unphased_geno[c1_indices], dtype=np.int8)
        geno2 = np.ascontiguousarray(unphased_geno[c2_indices], dtype=np.int8)
        
        # 4. Run C Function
        stats = Stats()
        if len(geno1) > 0 and len(geno2) > 0:
            _LD_LIB.compute_ld_stats(
                geno1, len(geno1),
                geno2, len(geno2),
                n_individuals,
                ctypes.byref(stats)
            )
            summary_results.append({
                'Time': t,
                'n_pairs': stats.n_pairs,
                'Value': stats.mean,
                'median_r2': stats.median,
                'sd_r2': stats.std,
                'min_r2': stats.min,
                'max_r2': stats.max
            })
        else:
            summary_results.append({
                'Time': t, 'n_pairs': 0, 
                'Value': 0, 'median_r2': 0, 'sd_r2': 0, 'min_r2': 0, 'max_r2': 0
            })
            
    return pd.DataFrame(summary_results)



##### Main Processing Function #####

def process_trees_file():
    for prefix in get_tree_prefixes():
        ts = tskit.load(f"{prefix}.trees")
        sample_id_dict = {t: ts.samples(time=t) for t in timepoints}
        
        # Oneway stats
        oneway_output(ts, prefix)
        
        # LD stats
        ld_summary = calculate_ld_stats(ts, sample_id_dict)
        ld_summary.to_csv(f"{prefix}_ld_summary.csv", index=False)

process_trees_file()
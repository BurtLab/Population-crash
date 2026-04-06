# Assessing alternative methods of using population genomic data to measure changes in population size

This directory contains the codes and data for producing the results in the article *Assessing alternative methods of using population genomic data to measure changes in population size*. All scripts were written and tested using **R 4.4.1**, **Python 3.12.3**, **Visual Studio Code 1.94.2**, and the **bash terminal** on **Ubuntu 24.04.1 LTS**.

## Project structure

- `code/`: Contains codes for simulations on HPC and Jupyter Notebook for analysis and plotting. 
- `data/`: Statistics outputs transfered from HPC that will be use for analysis in the Notebooks.
- `outputs/`: An empty directory for storage of images produced from the codes.

---

## Usage

### Running the scripts
- **hpc_*.sh**

All the shell scripts with the hpc_ prefix are for batch jobs on PBS HPC system. Transfer them along with all the Python and R scripts to the `$HOME` directory on HPC, then submit the jobs by `qsub` command. E.g. `qsub hpc_LD.sh`.

### Dependencies
You need to first set up conda on the remote HPC system. This step may vary on different HPC systems, please refer to instructions from your own institution.

In conda environment, install R by

```
conda install R
```

Then set up a virtual environment by
```
conda create -n msprime_env msprime tskit numpy pandas r-tidyverse r-reticulate 
```
This will create an environment `msprime_env` with all necessary modules and packages installed.

---

## Author
Laiyin Zhou
l.zhou24@imperial.ac.uk
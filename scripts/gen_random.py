import numpy as np

n_rows = snakemake.params.n_rows
n_cols = snakemake.params.n_cols
mean = snakemake.params.mean
sigma = snakemake.params.sigma
seed = snakemake.params.seed

np.random.seed(seed)
np.savetxt(snakemake.output, np.random.normal(mean, sigma, (n_rows, n_cols)), delimiter="\t")
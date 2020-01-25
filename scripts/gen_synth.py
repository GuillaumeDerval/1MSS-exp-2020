import json
import sys
import numpy as np
import os.path

def gen(output_path, file_prefix_name, n_to_generate, p_row_sol, p_col_sol,
        n_rows, n_cols, mean_high_distribution, mean_low_distribution,
        sigma_high_distribution, sigma_low_distribution):
    BASE_SEED = 828927

    np.random.seed(BASE_SEED)
    seeds = np.random.randint(0, 10000000, (n_to_generate,))

    base_description = {
        "p_row_sol": p_row_sol,
        "p_col_sol": p_col_sol,
        "n_rows": n_rows,
        "n_cols": n_cols,
        "mean_high_distribution": mean_high_distribution,
        "mean_low_distribution": mean_low_distribution,
        "sigma_low_distribution": sigma_low_distribution,
        "sigma_high_distribution": sigma_high_distribution
    }

    for file_entry in range(n_to_generate):
        seed = seeds[file_entry]
        np.random.seed(file_entry)
        selected_rows = [i for i, p in enumerate(np.random.rand(n_rows)) if p <= p_row_sol]
        selected_cols = [i for i, p in enumerate(np.random.rand(n_cols)) if p <= p_col_sol]

        matrix = np.random.normal(mean_low_distribution, sigma_low_distribution, (n_rows, n_cols))
        intern_matrix = np.random.normal(mean_high_distribution, sigma_high_distribution, (len(selected_rows), len(selected_cols)))
        matrix[np.ix_(selected_rows,selected_cols)] = intern_matrix

        filename_matrix = os.path.join(output_path, f"{file_prefix_name}{file_entry}.tsv")
        filename_description = os.path.join(output_path, f"{file_prefix_name}{file_entry}_description.json")
        np.savetxt(filename_matrix, matrix, fmt='%.4f', delimiter="\t")

        description = {"id": file_entry, "chosen_rows": selected_rows, "chosen_cols": selected_cols}
        description.update(base_description)
        json.dump(description, open(filename_description, 'w'))

    filename_description = os.path.join(output_path, f"{file_prefix_name}description.json")
    base_description["n_generated"] = n_to_generate
    json.dump(base_description, open(filename_description, 'w'))


if __name__ == '__main__':
    try:
        output_path = snakemake.params.output_path
        file_prefix_name = snakemake.params.file_prefix_name
        n_to_generate = snakemake.params.n_to_generate
        p_row_sol = snakemake.params.p_row_sol
        p_col_sol = snakemake.params.p_col_sol
        n_rows = snakemake.params.n_rows
        n_cols = snakemake.params.n_cols
        mean_high_distribution = snakemake.params.mean_high_distribution
        mean_low_distribution = snakemake.params.mean_low_distribution
        sigma_low_distribution = snakemake.params.sigma_low_distribution
        sigma_high_distribution = snakemake.params.sigma_high_distribution
    except:
        (output_path, file_prefix_name, n_to_generate, p_row_sol, p_col_sol,
         n_rows, n_cols, mean_high_distribution,
         mean_low_distribution, sigma_high_distribution,
         sigma_low_distribution) = sys.argv[1:]

    n_to_generate, n_rows, n_cols = [int(x) for x in [n_to_generate, n_rows, n_cols]]
    p_row_sol, p_col_sol, mean_high_distribution, mean_low_distribution, sigma_low_distribution, sigma_high_distribution = [float(x) for x in [p_row_sol, p_col_sol, mean_high_distribution, mean_low_distribution, sigma_low_distribution, sigma_high_distribution]]

    gen(output_path, file_prefix_name, n_to_generate, p_row_sol, p_col_sol,
        n_rows, n_cols, mean_high_distribution,
        mean_low_distribution, sigma_low_distribution,
        sigma_high_distribution)

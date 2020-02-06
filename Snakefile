configfile: "config.yaml"

#
# Produce everything
#
rule all:
    input:
        "results/synthetic_small.db"

#
# SYNTHETIC SMALL (COMPLETE)
#

def produce_synthetic_small_sql_gen_ids():
    out = []
    for size in config["synth_small_sizes"]:
        for mean, sigma in config["synth_small_fill_with"]:
            for methodDirect in config["ctr_methods"]:
                for branching in config["branch_methods"]:
                    for dataset_id in range(config["synth_small_n_instances_per_size"]):
                        out.append(f"{size}_{mean}_{sigma}_{dataset_id}_{methodDirect}_{branching}")
    return out

rule produce_synthetic_small_sql:
    input:
         expand("results/synthetic_small/{d}.txt", d=produce_synthetic_small_sql_gen_ids())
    output:
        "results/synthetic_small.db"
    params:
        format=r"(?:.*[^0-9])?(?P<size>[0-9]+)_(?P<mean>[\-0-9\.]+)_(?P<sigma>[0-9\.]+)_(?P<dataset_id>[0-9]+)_(?P<ctr>[^_]+)_(?P<branching>[^_]+)\.txt",
        params=[("size", int), ("mean", float), ("sigma", float), ("dataset_id", int), ("ctr", str), ('branching', str)]
    script:
        "scripts/run_to_sqlite.py"

rule start_synthetic_small:
    input:
        file="data/synthetic_small/{dataset}.tsv",
        jar=config["jarpath"]
    output: "results/synthetic_small/{dataset}_{ctr}_{branching}.txt"
    group: "synth_{dtype}_{dataset}"
    shell: "java -cp {input.jar} oscar1mss.runners.CompleteRunner {input.file} 7200 {wildcards.ctr} {wildcards.ctr} {wildcards.branching} > {output}"

rule produce_synthetic_small:
    output: "data/synthetic_small/{size}_{mean}_{sigma}_{dataset_id}.tsv"
    run:
        import numpy as np
        np.random.seed(7893397+int(wildcards.dataset_id)*int(wildcards.size))
        np.savetxt(output[0], np.random.normal(float(wildcards.mean), float(wildcards.sigma), (int(wildcards.size), int(wildcards.size))), delimiter="\t")

#
# SYNTHETIC BIG (LNS)
#

def produce_synthetic_big_sql_gen_ids():
    out = []
    for size, psol, timeout, first_n_instances_considered in config["synth_big_sizes"]:
        for method_direct in config["ctr_methods"]:
            for branching in config["branch_methods"]:
                for dataset in range(first_n_instances_considered):
                    out.append(f"{size}_{psol}_{dataset}_{method_direct}_{branching}_{timeout}")
    return out

rule produce_synthetic_big_sql:
    input:
         expand("results/synthetic_big/{d}.txt", d=produce_synthetic_big_sql_gen_ids())
    output:
        "results/synthetic_big.db"
    params:
        format=r"(?:.*[^0-9])?(?P<n_rows>[0-9]+)x(?P<n_cols>[0-9]+)_(?P<p_sol>[0-9]+)_(?P<dataset_id>[0-9]+)_(?P<ctr>[^_]+)_(?P<branching>[^_]+)_(?P<max_time>[0-9]+)\.txt",
        params=[("n_rows", int), ("n_cols", int), ("p_sol", int), ("dataset_id", int), ("ctr", str), ('branching', str), ('max_time', int)]
    script:
        "scripts/run_to_sqlite.py"

rule start_synthetic_big:
    input:
        file="data/synthetic_big/{dtype}/{dataset}.tsv",
        jar=config["jarpath"]
    output: "results/synthetic_big/{dtype}_{dataset}_{ctr}_{branching}_{timeout}.txt"
    group: "synth_{dtype}_{dataset}"
    shell: "java -cp {input.jar} oscar1mss.runners.CompleteRunner {input.file} {wildcards.timeout} {wildcards.ctr} {wildcards.ctr} {wildcards.branching} > {output}"

rule produce_synthetic_big:
    input: "scripts/gen_synth_big.py"
    output: expand("data/synthetic_big/{{n_rows}}x{{n_cols}}_{{p_sol}}/{i}.tsv", i=range(config["synth_big_n_instances_per_size"]))
    params:
        output_path = "data/synthetic_big/{n_rows}x{n_cols}_{p_sol}",
        file_prefix_name = "",
        n_to_generate = int(config["synth_big_n_instances_per_size"]),
        p_row_sol = lambda wildcards, output: float(wildcards.p_sol)/100.0,
        p_col_sol = lambda wildcards, output: float(wildcards.p_sol)/100.0,
        n_rows = lambda wildcards, output: int(wildcards.n_rows),
        n_cols = lambda wildcards, output: int(wildcards.n_cols),
        mean_high_distribution = 0.5,
        mean_low_distribution = -0.5,
        sigma_low_distribution = 1,
        sigma_high_distribution = 1
    script: "scripts/gen_synth_big.py"


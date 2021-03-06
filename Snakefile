configfile: "config.yaml"

#
# Produce everything
#
rule all:
    input:
        "results/synthetic_small.db",
        "results/synthetic_big.db",
        "results/lagrange_convergence.json"

#
# SYNTHETIC SMALL (COMPLETE)
#

def produce_synthetic_small_sql_gen_ids():
    out = []
    for size in config["synth_small_sizes"]:
        for mean, sigma in config["synth_small_fill_with"]:
            for method, branching in config["methods"]:
                if branching != "mip":
                    for dataset_id in range(config["synth_small_n_instances_per_size"]):
                        out.append(f"{size}_{mean}_{sigma}_{dataset_id}_{method}_{branching}")
    return out

rule produce_synthetic_small_sql:
    input:
         expand("results/synthetic_small/{d}.txt", d=produce_synthetic_small_sql_gen_ids())
    output:
        "results/synthetic_small.db"
    params:
        format=r"(?:.*[^0-9])?(?P<size>[0-9]+)_(?P<mean>[\-0-9\.]+)_(?P<sigma>[0-9\.]+)_(?P<dataset_id>[0-9]+)_(?P<method>[^_]+)_(?P<branching>[^_]+)\.txt",
        params=[("size", int), ("mean", float), ("sigma", float), ("dataset_id", int), ("method", str), ('branching', str)]
    script:
        "scripts/run_to_sqlite.py"


rule start_synthetic_small:
    input:
        file="data/synthetic_small/{dataset}.tsv",
        jar=ancient(config["jarpath"])
    output:
        final="results/synthetic_small/{dataset}_{method}_{branching}.txt",
        tmp=temp("tmp/results/synthetic_small/{dataset}_{method}_{branching}.txt")
    group: "synth_{dataset}"
    resources:
            mem_mb=2200
    shell:
        "java -Xmx2G -cp {input.jar} oscar1mss.runners.CompleteRunner {input.file} 72000 {wildcards.method} {wildcards.branching} > {output.tmp};"
        "cp {output.tmp} {output.final};"
        "echo 'FINAL' >> {output.final}"

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
    for size, psol, timeout, first_n_instances_considered, lns_timeout in config["synth_big_sizes"]:
        for method, branching in config["methods"]:
            for dataset in range(first_n_instances_considered):
                out.append(f"{size}_{psol}_{dataset}_{method}_{branching}_{timeout}_{lns_timeout}")
    return out

rule produce_synthetic_big_sql:
    input:
         expand("results/synthetic_big/{d}.txt", d=produce_synthetic_big_sql_gen_ids())
    output:
        "results/synthetic_big.db"
    params:
        format=r"(?:.*[^0-9])?(?P<n_rows>[0-9]+)x(?P<n_cols>[0-9]+)_(?P<p_sol>[0-9]+)_(?P<dataset_id>[0-9]+)_(?P<method>[^_]+)_(?P<branching>[^_]+)_(?P<max_time>[0-9]+)_(?P<lns_timeout>[0-9]+)\.txt",
        params=[("n_rows", int), ("n_cols", int), ("p_sol", int), ("dataset_id", int), ("method", str), ('branching', str), ('max_time', int), ('lns_timeout', int)]
    script:
        "scripts/run_to_sqlite.py"

rule start_synthetic_big:
    input:
        file="data/synthetic_big/{dtype}/{dataset}.tsv",
        jar=ancient(config["jarpath"])
    output:
        final="results/synthetic_big/{dtype}_{dataset}_{method}_{branching}_{timeout}_{lns_timeout}.txt",
        tmp=temp("tmp/results/synthetic_big/{dtype}_{dataset}_{method}_{branching}_{timeout}_{lns_timeout}.txt")
    group: "synth_{dtype}_{dataset}"
    wildcard_constraints:
        branching="(static|dynamic)"
    resources:
        mem_mb=5500
    shell:
        "java -Xmx5G -cp {input.jar} oscar1mss.runners.LNSRunner2 {input.file} {wildcards.timeout} {wildcards.method} {wildcards.branching} {wildcards.lns_timeout} > {output.tmp};"
        "cp {output.tmp} {output.final};"
        "echo 'FINAL' >> {output.final}"

rule start_synthetic_big_mip:
    input:
        file="data/synthetic_big/{dtype}/{dataset}.tsv",
        mipscript=ancient(config["mippath"])
    output:
        final="results/synthetic_big/{dtype}_{dataset}_mip-{method}_mip_{timeout}_{lns_timeout}.txt",
        tmp=temp("tmp/results/synthetic_big/{dtype}_{dataset}_mip-{method}_mip_{timeout}_{lns_timeout}.txt")
    group: "synth_{dtype}_{dataset}"
    resources:
        mem_mb=5500
    shell:
        "export OPENBLAS_NUM_THREADS=1;"
        "python {input.mipscript} {input.file} {wildcards.method} {wildcards.timeout} > {output.tmp};"
        "cp {output.tmp} {output.final};"
        "echo 'FINAL' >> {output.final}"

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
        mean_high_distribution = -0.01,
        mean_low_distribution = -0.1,
        sigma_low_distribution = 1,
        sigma_high_distribution = 1
    script: "scripts/gen_synth_big.py"

#
# Expes Lagrange convergence
#

def run_lagrange_ids():
    all = set()
    for size in config["lagrange_sizes"]:
        for mean, sigma in config["lagrange_fill_with"]:
            for dataset_id in range(config["lagrange_n_instances_per_size"]):
                all.add(f"{size}_{mean}_{sigma}_{dataset_id}")
    return list(all)

rule sumup_lagrange:
    input: expand("results/lagrange_convergence/{d}.json", d=run_lagrange_ids())
    output: "results/lagrange_convergence.json"
    run:
        import json
        out = {}
        for f in input:
            with open(f, 'r') as ff:
                inp = json.load(ff)
                out[f] = inp
        json.dump(out, open(str(output), 'w'))

rule run_lagrange_convergence:
    input: "data/synthetic_small/{size}_{mean}_{sigma}_{dataset_id}.tsv"
    output: "results/lagrange_convergence/{size}_{mean}_{sigma}_{dataset_id}.json"
    script: "scripts/exp_lagrange_convergence.py"

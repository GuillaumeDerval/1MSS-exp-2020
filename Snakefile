configfile: "config.yaml"

#
# Produce everything
#
rule all:
    input:
        "results/find_shortest_among_all_branching_ctr_possibilities.db"

#
# Convert runs to SQLite
#
def find_shortest_among_all_branching_ctr_possibilities_gen_ids():
    out = []
    for size, psol, timeout, first_n_instances_considered in config["synth_test_complete"]:
        for methodDirect in config["ctr_methods"]:
            for methodTranspose in config["ctr_methods"]:
                for branching in config["branch_methods"]:
                    for dataset in range(first_n_instances_considered):
                        out.append("{size}_{psol}_{dataset}_{methodDirect}_{methodTranspose}_{branching}_{timeout}".format(
                            size=size, psol=psol, dataset=dataset,
                            methodDirect=methodDirect, methodTranspose=methodTranspose,
                            branching=branching, timeout=timeout
                        ))
    return out

rule find_shortest_among_all_branching_ctr_possibilities:
    input:
         expand("results/synthetic/{d}.txt", d=find_shortest_among_all_branching_ctr_possibilities_gen_ids())
    output:
        "results/find_shortest_among_all_branching_ctr_possibilities.db"
    script:
        "scripts/run_to_sqlite.py"

#
# Run things
#

rule produce_synthetic:
    input: "scripts/gen_synth.py"
    output: expand("data/{{n_rows}}x{{n_cols}}_{{p_sol}}/{i}.tsv", i=range(config["n_instances_per_param"]))
    params:
        output_path = "data/{n_rows}x{n_cols}_{p_sol}",
        file_prefix_name = "",
        n_to_generate = int(config["n_instances_per_param"]),
        p_row_sol = lambda wildcards, output: float(wildcards.p_sol)/100.0,
        p_col_sol = lambda wildcards, output: float(wildcards.p_sol)/100.0,
        n_rows = lambda wildcards, output: int(wildcards.n_rows),
        n_cols = lambda wildcards, output: int(wildcards.n_cols),
        mean_high_distribution = 0.5,
        mean_low_distribution = -0.5,
        sigma_low_distribution = 1,
        sigma_high_distribution = 1
    script:
        "scripts/gen_synth.py"

rule start_synthetic:
    input:
        file="data/{dtype}/{dataset}.tsv",
        jar=config["jarpath"]
    output: "results/synthetic/{dtype}_{dataset}_{methodDirect}_{methodTranspose}_{branching}_{timeout}.txt"
    group: "synth_{dtype}_{dataset}"
    shell: "java -cp {input.jar} oscar1mss.runners.CompleteRunner {input.file} {wildcards.timeout} {wildcards.methodDirect} {wildcards.methodTranspose} {wildcards.branching} > {output}"
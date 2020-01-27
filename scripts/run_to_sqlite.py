import sqlite3
import os
import sys
sys.path.insert(0, os.path.join(snakemake.scriptdir, ".."))
print(snakemake.scriptdir)
from scripts.extract_data_run import extract_params_from_name, extract_sols, extract_runs


def gen_db(dbpath, runs_path):
    try:
        os.unlink(dbpath)
    except:
        pass

    conn = sqlite3.connect(dbpath)
    c = conn.cursor()

    # Create table
    c.execute('''CREATE TABLE runs (n_rows int, n_cols int, p_sol int, dataset_id int,
                                    ctr_direct text, ctr_transpose text, branching text, max_time int, 
                                    runtime int, n_nodes int, n_fails int, completed boolean, n_sols int, 
                                    best_sol_time int, best_sol_node int, best_sol_obj real)''')
    c.execute('''CREATE TABLE run_sols (run_rowid int, obj real, node int, time int)''')

    runs_add = []
    for f in runs_path:
        params = extract_params_from_name(f)
        runs = extract_runs(f)
        sols = extract_sols(f)

        runtime = 0
        n_nodes = 0
        n_fails = 0
        completed = False
        n_sols = len(sols)
        best_sol_time = 0
        best_sol_node = 0
        best_sol_obj = 0.0

        if len(runs) == 1:
            runtime = runs[0]["runtime"]
            n_nodes = runs[0]["n_nodes"]
            n_fails = runs[0]["n_fails"]
            completed = runs[0]["completed"]
            assert n_sols == runs[0]["n_sols"]

        if len(sols) != 0:
            best_sol_obj = sols[-1][0]
            best_sol_node = sols[-1][1]
            best_sol_time = sols[-1][2]

        entry = (params["n_rows"], params["n_cols"], params["p_sol"], params["dataset_id"],
                 params["ctr_direct"], params["ctr_transpose"], params["branching"], params["max_time"],
                 runtime, n_nodes, n_fails, completed, n_sols, best_sol_time, best_sol_node, best_sol_obj)

        c.execute('INSERT INTO runs (n_rows, n_cols, p_sol, dataset_id, ctr_direct, ctr_transpose, branching, max_time, runtime, n_nodes, n_fails, completed, n_sols, best_sol_time, best_sol_node, best_sol_obj) '
                  'VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)', entry)

        row_id = c.lastrowid
        c.executemany('INSERT INTO run_sols (run_rowid, obj, node, time) '
                      'VALUES (?,?,?,?)', [(row_id, obj, node, time) for obj, node, time in sols])


    conn.commit()
    conn.close()

gen_db(snakemake.output[0], snakemake.input)
#gen_db("test.db", ["../results/synthetic/20x20_30_0_lpbound-fast-no-recompute-filter-row_lpbound-fast-no-recompute-filter-row_dynamic_30.txt"])
import sqlite3
import os
import sys
sys.path.insert(0, os.path.join(snakemake.scriptdir, ".."))
from scripts.extract_data_run import extract_params_from_name, extract_sols, extract_runs


def gen_db(dbpath, runs_path, format, params_names_and_type):
    try:
        os.unlink(dbpath)
    except:
        pass

    conn = sqlite3.connect(dbpath)
    c = conn.cursor()

    type_to_sqltype = {
        int: "int",
        str: "text",
        float: "real"
    }
    columns_with_type = ", ".join("{} {}".format(x[0], type_to_sqltype[x[1]]) for x in params_names_and_type)
    columns_name = ", ".join("{}".format(x[0]) for x in params_names_and_type)

    # Create table
    c.execute(f'''CREATE TABLE runs({columns_with_type}, 
                                    runtime int, n_nodes int, n_fails int, completed boolean, n_sols int, 
                                    best_sol_time int, best_sol_node int, best_sol_obj real)''')
    c.execute('''CREATE TABLE run_sols (run_rowid int, obj real, node int, time int)''')

    runs_add = []
    for f in runs_path:
        params = extract_params_from_name(f, format, params_names_and_type)
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

        entry = [params[x] for x, _ in params_names_and_type] + [runtime, n_nodes, n_fails, completed, n_sols, best_sol_time, best_sol_node, best_sol_obj]

        c.execute(f'INSERT INTO runs ({columns_name}, runtime, n_nodes, n_fails, completed, n_sols, best_sol_time, best_sol_node, best_sol_obj) '
                  'VALUES ('+','.join(["?"] * (len(params_names_and_type) + 8))+')', entry)

        row_id = c.lastrowid
        c.executemany('INSERT INTO run_sols (run_rowid, obj, node, time) '
                      'VALUES (?,?,?,?)', [(row_id, obj, node, time) for obj, node, time in sols])


    conn.commit()
    conn.close()

gen_db(snakemake.output[0], snakemake.input, snakemake.params.format, snakemake.params.params)
import re

regex_extract_sol = re.compile(r"Best: ([0-9.]+) / [0-9.]+ ----- ([0-9]+) ----- ([0-9.]+)")
regex_extract_run = re.compile(r"""nNodes: ([0-9]+)
nFails: ([0-9]+)
time\(ms\): ([0-9]+)
completed: (true|false)
timeInTrail: ([0-9]+)
nSols: ([0-9]+)""", re.MULTILINE)


def extract_params_from_name(fname,
                             format=r"(?:.*[^0-9])?(?P<n_rows>[0-9]+)x(?P<n_cols>[0-9]+)_(?P<p_sol>[0-9]+)_(?P<dataset_id>[0-9]+)_(?P<ctr_direct>[^_]+)_(?P<ctr_transpose>[^_]+)_(?P<branching>[^_]+)_(?P<max_time>[0-9]+)\.txt",
                             params=(("n_rows", int), ("n_cols", int), ("p_sol", int), ("dataset_id", int), ("ctr_direct", str), ('ctr_transpose', str), ('branching', str), ('max_time', str))):
    regex_extract_params = re.compile(format)
    result = regex_extract_params.match(fname).groupdict()
    out = {x: t(result[x]) for x,t in params}
    return out

def extract_sols(fname):
    data = open(fname).read().strip().split("\n")
    bests = []
    for line in data:
        match = regex_extract_sol.match(line)
        if match:
            bests.append((float(match.group(1)), int(match.group(2)), int(float(match.group(3))*1000.0)))
    return bests

def extract_runs(fname):
    data = open(fname).read().strip()
    out = []
    for x in regex_extract_run.findall(data):
        out.append({
            "n_nodes": int(x[0]),
            "n_fails": int(x[1]),
            "runtime": int(x[2]),
            "completed": x[3] == "true",
            "timeInTrail": int(x[4]),
            "n_sols": int(x[5])
        })
    return out
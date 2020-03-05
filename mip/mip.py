import argparse
from gurobipy import *
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('file', type=argparse.FileType('r'))
parser.add_argument('type', type=str, choices=["strong", "bigm"])
parser.add_argument('timeout_in_seconds', type=int)
args = parser.parse_args()

inputfile = args.file
type = args.type
timeout = args.timeout_in_seconds


data = [[float(y) for y in x.split("\t")] for x in inputfile.read().split("\n") if x != ""]

nRows = len(data)
nCols = len(data[0])

if nRows >= 2000 and type == "strong":
    print("This will consume too much memory. Stopping.")
    exit(0)

rows_idx = range(0, nRows)
cols_idx = range(0, nCols)

cost = {(i,j): -data[i][j] for i in rows_idx for j in cols_idx}

m = Model('mymip')
m.setParam(GRB.Param.OutputFlag, False)
m.setParam(GRB.Param.Threads, 1)
m.setParam(GRB.Param.TimeLimit, timeout)

if type == "strong":
    cells = m.addVars(rows_idx, cols_idx, name="s", obj=cost, lb=0.0, ub=1.0, vtype=GRB.BINARY)
    rows = m.addVars(rows_idx, name="r", lb=0.0, ub=1.0, vtype=GRB.BINARY)
    cols = m.addVars(cols_idx, name="c", lb=0.0, ub=1.0, vtype=GRB.BINARY)

    m.addConstrs((rows[i] + cols[j] <= cells[i,j] + 1 for i in rows_idx for j in cols_idx))
    m.addConstrs((cells[i,j] <= rows[i] for i in rows_idx for j in cols_idx))
    m.addConstrs((cells[i,j] <= cols[j] for i in rows_idx for j in cols_idx))
elif type == "bigm":
    row_used = m.addVars(rows_idx, name="r", lb=0.0, ub=1.0, vtype=GRB.BINARY)
    col_used = m.addVars(cols_idx, name="c", lb=0.0, ub=1.0, vtype=GRB.BINARY)
    row_cont = m.addVars(rows_idx, name="v", lb=0.0, obj=-1.0)

    def sum_row(i):
        return quicksum([col_used[j] * data[i][j] for j in cols_idx])

    def bigMPlus(i):
        s = 0.0
        for j in cols_idx:
            s += max(0, data[i][j])
        return s

    def bigMMinus(i):
        s = 0.0
        for j in cols_idx:
            s -= min(0, data[i][j])
        return s

    m.addConstrs((row_cont[i] <= row_used[i] * bigMPlus(i) for i in rows_idx), "a")
    m.addConstrs((row_cont[i] <= sum_row(i) + ((1.0 - row_used[i]) * bigMMinus(i)) for i in rows_idx), "b")
else:
    print("Unknown model type")
    exit(1)


def mycallback(model, where):
  if where == GRB.Callback.MIPSOL:
    best = -model.cbGet(GRB.Callback.MIPSOL_OBJBST)
    node = int(model.cbGet(GRB.Callback.MIPSOL_NODCNT))
    time = model.cbGet(GRB.Callback.RUNTIME)
    print(f"Best: {best} / 0.0 ----- {node} ----- {time}")

m.optimize(mycallback)

time = m.getAttr(GRB.Attr.Runtime)
nodes = int(m.getAttr(GRB.Attr.NodeCount))
completed = "true" if m.getAttr(GRB.Attr.ObjBound) == m.getAttr(GRB.Attr.ObjVal) else "false"
nsols = m.getAttr(GRB.Attr.SolCount)
print(f"""nNodes: {nodes}
nFails: 0
time(ms): {time}
completed: {completed}
timeInTrail: 0
nSols: {nsols}""")

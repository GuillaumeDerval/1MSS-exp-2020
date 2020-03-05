import argparse
import json

from gurobipy import *
import numpy as np

try:
    inputfile = open(str(snakemake.input), "r")
    outputfile = open(str(snakemake.output), "w")
except:
    parser = argparse.ArgumentParser()
    parser.add_argument('input', type=argparse.FileType('r'))
    parser.add_argument('output', type=argparse.FileType('w'))
    args = parser.parse_args()

    inputfile = args.input
    outputfile = args.output


data = [[float(y) for y in x.split("\t")] for x in inputfile.read().split("\n") if x != ""]

def find_optimum(data):
    nRows = len(data)
    nCols = len(data[0])

    rows_idx = range(0, nRows)
    cols_idx = range(0, nCols)

    cost = {(i, j): -data[i][j] for i in rows_idx for j in cols_idx}

    m = Model('mymip')
    #m.setParam(GRB.Param.OutputFlag, False)
    m.setParam(GRB.Param.Threads, 1)

    cells = m.addVars(rows_idx, cols_idx, name="s", obj=cost, lb=0.0, ub=1.0, vtype=GRB.CONTINUOUS)
    rows = m.addVars(rows_idx, name="r", lb=0.0, ub=1.0, vtype=GRB.CONTINUOUS)
    cols = m.addVars(cols_idx, name="c", lb=0.0, ub=1.0, vtype=GRB.CONTINUOUS)

    m.addConstrs((rows[i] + cols[j] <= cells[i,j] + 1 for i in rows_idx for j in cols_idx))
    m.addConstrs((cells[i,j] <= rows[i] for i in rows_idx for j in cols_idx))
    m.addConstrs((cells[i,j] <= cols[j] for i in rows_idx for j in cols_idx))
    m.optimize()

    return -m.objVal

def run_lagrange(data):
    #optimum = find_optimum(data)
    data = np.array(data)
    nRows = len(data)
    nCols = len(data[0])
    alpha = np.random.uniform(0, 1, data.shape)
    gamma = np.random.uniform(0, 1, data.shape)

    for i in range(nRows):
        for j in range(nCols):
            if data[i,j] > 0:
                gamma[i,j] = 0
            else:
                alpha[i,j] = 0

    def compute_r_c_x():
        ub = 0.0
        r = np.zeros((nRows, ))
        c = np.zeros((nCols, ))
        x = np.zeros((nRows, nCols))

        for i in range(nRows):
            contribution = 0.0
            for j in range(nCols):
                ub += gamma[i,j]
                contribution += alpha[i,j] - gamma[i,j]
            if contribution > 0:
                r[i] = 1
                ub += contribution

        for j in range(nCols):
            contribution = 0.0

            for i in range(nRows):
                contribution -= gamma[i,j]
                cell_contribution = data[i,j] - alpha[i,j] + gamma[i,j]
                if cell_contribution > 0:
                    contribution += cell_contribution
                    x[i,j] = 1

            if contribution > 0:
                ub += contribution
                c[j] = 1
            else:
                for i in range(nRows):
                    x[i, j] = 0

        return ub, r, c, x

    def iteration(mu):
        ub, r, c, x = compute_r_c_x()
        for i in range(nRows):
            for j in range(nCols):
                if data[i, j] > 0:
                    alpha[i, j] = max(0.0, alpha[i, j]-mu*(r[i]-x[i, j]))
                else:
                    gamma[i, j] = max(0.0, gamma[i, j]-mu*(x[i, j]+1-r[i]-c[j]))
        return ub

    mu = 1.0
    out = []
    for itr in range(500):
        ub = iteration(mu)
        out.append(ub)
        mu*=0.95

    return out

out = run_lagrange(data)
json.dump(out, outputfile)
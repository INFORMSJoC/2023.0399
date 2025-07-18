import pandas as pd
import numpy as np
from tabulate import tabulate
from natsort import os_sorted
import sys; sys.path.insert(0, '..')
from runAll import *

projectDir = "../"

instanceRoot = projectDir + "data/CSP"

instancesSuffixes = ["NewHard-20to40-400-BPP", "NewHard-20to40-200-BPP", "NewHard-20to40-600-BPP", "AI601", "AI802", "Hard", "ANI201", "ANI402", "ANI600",  "ANI801", "ANI1002", "Random", "FalkenauerT", "FalkenauerU", "Schwerin", "Waescher", "AI202", "AI403", "AI1003", "Scholl", "IrnichAA", "IrnichAB", "IrnichBA", "IrnichBB"]

checkerCSP = projectDir + "src/checkerCSP"

alg = projectDir + "src/CSPSolver"
logs = [projectDir + "src/CSPSolver/out/" + suff + ".log" for suff in instancesSuffixes]
features = ["RF", "CRF", "Waste", "Cost", "Cut", "Historic", "Splay", "Ineq"] #, Best

for suf in instancesSuffixes:
    instanceDir = instanceRoot + "/" + suf
    runAll(alg, instanceDir, checkerCSP, False, features, softClean=False)

tables = []
for log in logs:
    table = pd.read_csv(log, sep=',', usecols=[i for i in range(7)] + [10, 11, 12, 13, 15, 16, 17, 18],header=None, names = ['Name', 'Z', 'RootObj', 'Time', 'Optimal', 'genCols', 'genCuts', 'IntR', 'FratR', 'PTime', 'LPTime', 'Where', 'Branches', 'PCall1', 'PCall2'])
    tables.append(table)

res = []
def dec2(a):
    if a < 0.01:
        return round(a * 1000)/1000
    if a < 0.1:
        return round(a * 100)/100
    return round(a * 10)/10

for h in range(len(tables)):
    table = tables[h]
    table.loc[:, 'OptR'] = (table.Where == 0) | (table.Where == 1)
    table.loc[:, 'OptH'] = (table.Where == 2) & (table.Branches <= 0)
    table.loc[:, 'Opt5'] = (table.Where == 2) & (table.Branches <= 5) & (table.Branches > 0)
    table.loc[:, 'OptM'] = (table.Where == 2) & (table.Branches > 5)
    means = table.mean()
    means.IntR *= 100
    means = means.apply(dec2)
    res.append([instancesSuffixes[h], table.shape[0], table.Optimal.eq('OPT').sum(), table.OptM.sum(), means.Time, means.PTime, means.LPTime, means.genCols, means.genCuts, means.IntR, means.PCall1, means.PCall2, dec2(table.genCols.sum() / table.PCall2.sum())])

res = os_sorted(res, key=lambda x: x[0])
res = [["Name", "Total", "Opt",  "OptM", "Time",  "PTime", "LPTime", "Cols", "Cuts", "IntR", "PCall1", "PCall2", "PCall3"]] + res
#print(tabulate(res))

res = np.array([np.array(xi) for xi in res])
df = pd.DataFrame(res[1:, 1:], index=res[1:,0], columns=res[0, 1:])
print(df.to_latex())

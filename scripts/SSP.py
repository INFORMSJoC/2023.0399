import pandas as pd
import numpy as np
from natsort import os_sorted
from tabulate import tabulate
import sys; sys.path.insert(0, '..')
from runAll import *

projectDir = "../"

instanceRoot = projectDir + "data/SSP/"

instancesSuffixes= ["A1", "A2", "A3", "B", "AI202", "ANI201", "AI403", "ANI402", "GI125", "GI250", "GI500", "GI750", "GI1000", "Hard", "FalkenauerT", "FalkenauerU", "Schwerin", "Waescher", "Scholl", "ANI600", "ANI801", "AI601", "AI802", "AI1003", "ANI1002"] 

checkerBPP = projectDir + "src/checkerSSP"

alg = projectDir + "src/SSPSolver"
logs = [projectDir + "src/SSPSolver/out/" + suff + ".log" for suff in instancesSuffixes]
features = ["RF", "CRF", "Waste", "Cost", "Cut", "Historic", "Splay", "Ineq"]

for suf in instancesSuffixes:
    instanceDir = instanceRoot + "/" + suf
    runAll(alg, instanceDir, checkerBPP, False, features, softClean=False)

tables = []
for log in logs:
    table = pd.read_csv(log, sep=',', usecols=[i for i in range(7)] + [11, 12, 13, 14],header=None, names = ['Name', 'Z', 'RootObj', 'Time', 'Optimal', 'genCols', 'genCuts', 'IntR', 'FratR', 'PTime', 'LPTime'])
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
    means = table.mean()
    means.IntR *= 100
    means = means.apply(dec2)
    res.append([instancesSuffixes[h], table.shape[0], table.Optimal.eq('OPT').sum(), means.Time, means.genCols, means.genCuts, means.IntR, means.FratR, means.PTime, means.LPTime])

res = os_sorted(res, key=lambda x: x[0])
res = [["Name", "Total", "Opt", "Time", "Cols", "Cuts", "IntR", "FratR", "PTime", "LPTime"]] + res
#print(tabulate(res))

res = np.array([np.array(xi) for xi in res])
df = pd.DataFrame(res[1:, 1:], index=res[1:,0], columns=res[0, 1:])
df = df.reset_index().rename(columns= {'index' : 'Name'})[['Name', 'Total', 'Opt', 'Time', 'Cols']]
df2 = pd.read_csv("SSP.csv", sep=' ')
df = df.astype('str').merge(df2.astype('str'), on=['Name', 'Total'], how='left')[['Name', 'Total', 'OptC', 'TimeC', 'OptB', 'TimeB', 'OptA', 'TimeA', 'Opt', 'Time', 'Cols']]
for i in range(df.shape[0]):
    colsO = ['OptC', 'OptA', 'OptB', 'Opt']
    mx = str(df.loc[i, colsO].replace({'–':'0'}).astype('int').max())
    for j in colsO:
        if df.loc[i, j] == mx:
            df.loc[i, j] = "\textbf{" + mx + "}"

    colsT = ['TimeC', 'TimeA', 'TimeB', 'Time']
    mx = str(df.loc[i, colsT].replace({'–':'inf'}).astype('float').min())
    for j in colsT:
        if df.loc[i, j] != '–' and df.loc[i, j] == mx:
            df.loc[i, j] = "\textbf{" + mx + "}"

df = df.set_index('Name')
print(df.to_latex(escape=False))

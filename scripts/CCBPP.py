import pandas as pd
import numpy as np
from math import ceil
from natsort import os_sorted
import sys; sys.path.insert(0, '..')
from runAll import *

projectDir = "../"

instanceRoot = projectDir + "data/CCBPP/"

instancesSuffixes = ["Borges"]

checkerBPP = projectDir + "src/checkerCCBPP"

alg = projectDir + "src/CCBPPSolver"
logs = [projectDir + "src/CCBPPSolver/out/" + suff + ".log" for suff in instancesSuffixes]
features = ["RF", "CRF", "Waste", "Cost", "Cut", "Historic", "Splay"]

for suf in instancesSuffixes:
    instanceDir = instanceRoot + "/" + suf
    runAll(alg, instanceDir, checkerBPP, False, features, softClean=False)

tables = []
for log in logs:
    table = pd.read_csv(log, sep=',', usecols=[i for i in range(7)],header=None, names = ['Name', 'Z', 'RootObj', 'Time', 'Optimal', 'genCols', 'genCuts'])
    tables.append(table)

names = ['Name', 'N', 'Q', 'C', 'W', 'RootObjY', 'ZY', 'genColsY', 'TimeY']
yulle = pd.read_csv("test1-short.out", sep = "\t", usecols = [0, 1, 2, 3, 4, 8, 9, 11, 21], names=names)
yulle.loc[:, 'Name'] = [i[:-1] for i in yulle.Name]
yulle.loc[:, 'OptimalY'] = np.where(yulle.TimeY < 900000, "OPT", "TLE")
yulle.loc[:, 'TimeY'] = yulle.TimeY / 1000

table = table.merge(yulle, on=['Name'])
table.loc[:, 'value'] = table.Name.apply(lambda x : int(x[1])*4 + int(x[3]))
gb = table.groupby('value')
tables = [gb.get_group(x).reset_index(drop=True) for x in gb.groups]
mask = {5 : 'Q10C2', 6 : 'Q10C3', 7 : 'Q10CX', 9 : 'Q25C2', 10 : 'Q25C3', 11 : 'Q25CX', 13 : 'Q50C2', 14 : 'Q50C3', 15 : 'Q50CX'}
instancesSuffixes = [mask[tables[i].iloc[0].value] for i in range(9)]

res = []
def dec2(a):
    if a < 0.01:
        return round(a * 1000)/1000
    if a < 0.1:
        return round(a * 100)/100
    return round(a * 10)/10

for h in range(len(tables)):
    table = tables[h]
    N = len(table)
    means = table.mean()
    means = means.apply(dec2)
    acumANI = 0
    acumANI2 = 0
    acumANIgap = 0
    for i in range(len(table)):
        rootObj = ceil(float(table["RootObj"][i]) - 1e-10)
        ZSol = int(table["Z"][i])
        acumANI += ZSol > rootObj
        acumANI2 += ZSol > rootObj + 1
        if ZSol > rootObj:
            acumANIgap += ZSol - float(table["RootObj"][i])
    res.append([instancesSuffixes[h], N, table.OptimalY.eq('OPT').sum(), means.TimeY, means.genColsY, table.Optimal.eq('OPT').sum(), means.Time, means.genCols, means.genCuts, acumANI, acumANI2, acumANIgap])

res = os_sorted(res, key=lambda x: x[0])
res = [["Name", "Total", "OptY", "TimeY", "ColsY", "Opt", "Time", "Cols", "Cuts", "ANI", "ANI2", "ANIgap"]] + res

res = np.array([np.array(xi) for xi in res])
df = pd.DataFrame(res[1:, 1:], index=res[1:,0], columns=res[0, 1:])
df = df.reset_index().astype('str').rename(columns= {'index' : 'Name'})[['Name', 'Total', 'ANI', 'OptY', 'TimeY', 'ColsY', 'Opt', 'Time', 'Cols', 'Cuts']]

for i in range(df.shape[0]):
    colsO = ['OptY', 'Opt']
    mx = str(df.loc[i, colsO].replace({'–':'0'}).astype('int').max())
    for j in colsO:
        if df.loc[i, j] == mx:
            df.loc[i, j] = "\textbf{" + mx + "}"

    colsT = ['TimeY', 'Time']
    mx = str(df.loc[i, colsT].replace({'–':'inf'}).astype('float').min())
    for j in colsT:
        if df.loc[i, j] != '–' and df.loc[i, j] == mx:
            df.loc[i, j] = "\textbf{" + mx + "}"

df.loc[:, 'C'] = [n[-1] for n in df.Name]
df.loc[:, 'Q'] = [n[1:3] for n in df.Name]
df = df[['Q', 'C', 'Total', 'ANI', 'OptY', 'TimeY', 'ColsY', 'Opt', 'Time', 'Cols', 'Cuts']]
df = df.set_index(['Q', 'C'])
print(df.to_latex(escape=False))

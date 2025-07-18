import pandas as pd
import numpy as np
from natsort import os_sorted
from tabulate import tabulate
import sys; sys.path.insert(0, '..')
from runAll import *

projectDir = "../"

instanceRoot = projectDir + "data/IPMS"

suff1 = ["2", "2.25", "2.5", "2.75", "3"]
instancesSuffixes = []
for s1 in suff1:
    for i in range(1, 8):
        instancesSuffixes.append("Ratio" + s1 + "/Ratio" + s1 + "-Class" + str(i))
     
checkerIPMS = projectDir + "src/checkerIPMS"

alg = projectDir + "src/IPMSSolver"
logs = [projectDir + "src/IPMSSolver/out/" + suff[suff.rfind('/') + 1: ] + ".log" for suff in instancesSuffixes]
features = ["RF", "CRF", "Waste", "Cost", "Cut", "Historic", "Splay"]

for suf in instancesSuffixes:
    instanceDir = instanceRoot + "/" + suf
    print(instanceDir)
    runAll(alg, instanceDir, checkerIPMS, False, features, softClean=False)


tables = []
for log in logs:
    table = pd.read_csv(log, sep=',', usecols=[i for i in range(10)],header=None, names = ['Name', 'W', 'Machines', 'Time', 'Optimal', 'Diff', 'Type', 'genCols', 'genCuts', 'iter'])
    tables.append(table)

res = []
def dec2(a):
    if a < 0.01:
        return round(a * 1000)/1000
    if a < 0.1:
        return round(a * 100)/100
    return round(a * 10)/10

tables2 = {}
for h in range(len(tables)):
    table = tables[h]
    name = instancesSuffixes[h]
    pref, cl = name.split('-')
    if cl[-1] == '6' or cl[-1] == '7':
        cl = "Class6-7"
    else:
        cl = "Class1-5"
    
    pref = pref + "-" + cl
    
    tab = table.copy()
    if pref in tables2:
        tab = pd.concat([tab, tables2[pref]], ignore_index=True)
    
    tables2[pref] = tab

for key in tables2:
    table = tables2[key]
    N = len(table)
    acumTime = 0
    acumOpt = 0
    acumCols = 0
    acumCuts = 0
    for i in range(len(table)):
        acumTime += float(table["Time"][i])
        acumOpt += table["Optimal"][i] == "OPT"
        acumCols += int(table["genCols"][i])
        acumCuts += int(table["genCuts"][i])
    
    res.append([key, N, acumOpt, dec2(acumTime / N), dec2(acumCols / N), dec2(acumCuts / N)])

res = os_sorted(res, key=lambda x: x[0])
res = [["Name", "Total", "Opt", "Time", "Cols", "Cuts"]] + res
print(tabulate(res))
res = np.array([np.array(xi) for xi in res])
df = pd.DataFrame(res[1:, 1:], index=res[1:,0], columns=res[0, 1:])
df.to_csv("IPMS.csv")

res = []
for h in range(len(tables)):
    table = tables[h]
    N = len(table)
    acumTime = 0
    acumOpt = 0
    acumCols = 0
    acumCuts = 0
    for i in range(len(table)):
        acumTime += float(table["Time"][i])
        acumOpt += table["Optimal"][i] == "OPT"
        acumCols += int(table["genCols"][i])
        acumCuts += int(table["genCuts"][i])
        
    res.append([instancesSuffixes[h], N, acumOpt, dec2(acumTime / N), dec2(acumCols / N), dec2(acumCuts / N)])

res = os_sorted(res, key=lambda x: x[0])
res = [["Name", "Total", "Opt", "Time", "Cols", "Cuts"]] + res
print(tabulate(res))


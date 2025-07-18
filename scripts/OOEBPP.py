import pandas as pd
from tabulate import tabulate
from natsort import os_sorted
import sys; sys.path.insert(0, '..')
from runAll import *

projectDir = "../"

instanceRoot = projectDir + "data/OOEBPP"

instancesSuffixes = ["R50", "R100", "R200", "R300", "R400", "R500", "R750", "R1000", "BENG", "CGCUT", "GCUT1-4","GCUT5-13","NGCUT", "CLASS", "HT"]

checkerBPP = projectDir + "src/checkerOOEBPP"

alg = projectDir + "src/OOEBPPSolver"
logs = [projectDir + "src/OOEBPPSolver/out/" + suff + ".log" for suff in instancesSuffixes]
features = ["RF", "CRF", "Waste", "Cost", "Cut", "Historic", "Splay"]

for suf in instancesSuffixes:
    instanceDir = instanceRoot + "/" + suf
    runAll(alg, instanceDir, checkerBPP, False, features, softClean=False)
    
tables = []
for log in logs:
    table = pd.read_csv(log, sep=',', usecols=[i for i in range(7)],header=None, names = ['Name', 'Z', 'RootObj', 'Time', 'Optimal', 'genCols', 'genCuts'])
    tables.append(table)

res = []
def dec2(a):
    if a < 0.001:
        return round(a * 10000)/10000
    if a < 0.01:
        return round(a * 1000)/1000
    if a < 0.1:
        return round(a * 100)/100
    return round(a * 10)/10

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

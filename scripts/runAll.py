import os
from glob import glob

def runAll(dirHeuristic, instancesDir, checker, clearLog, features = [], softClean = True, seed = -1):
    heuristicName = dirHeuristic[dirHeuristic.rfind('/') + 1:]
    heuristic = dirHeuristic + "/" + heuristicName
    
    checkerName = checker[checker.rfind("/checker") + 1:]
    os.system("(cd " + checker[:checker.rfind("/")] + "; make " +  checkerName + ")")
    
    os.system("cd " + dirHeuristic + "; make")

    if not os.path.exists(instancesDir):
        print("Not found the dir of instance")
        exit(1)

    if(clearLog):
        os.system("rm -r " + dirHeuristic + "/out")
    os.system("mkdir " + dirHeuristic + "/out")

    instancesDirName = instancesDir[instancesDir.rfind('/') + 1:]
    logFile = dirHeuristic + "/out/" + instancesDirName + ".log"
    logAll = open(logFile, "w")
    
    if len(features) > 0:
        lfeatures = " -" + " -".join(features)
        print(lfeatures)
    else:
        lfeatures = ""
    if seed != -1:
        lfeatures += " -seed " + str(seed)
    print(lfeatures)

    for instance in sorted(glob(os.path.join(instancesDir, "*"))):
        print ("Running instance", instance)
        instanceName = instance[instance.rfind('/') + 1:]
        instanceLog = dirHeuristic + "/out/" + instanceName + ".log"
        if not os.path.exists(instanceLog):
            os.system(heuristic + " " + instance + " " + checker + lfeatures)
        logFileCur = open(instanceLog, "r")
        logAll.write(instanceName + ',' + logFileCur.readline())
        logAll.flush()
        logFileCur.close()
    
    if softClean:
        os.system("(cd " + dirHeuristic + " ; make softClean)")
    return logFile

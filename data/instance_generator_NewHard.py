#! /usr/bin/python3

'''
A generator of hard instances for the BPP.

The instances generate are made cutting the bin, that is,
the value of the optimal solution is equal to lower bound L1

Info:
 - The items are in range [0.25, 0.5]
 - I do not know if the algorithm works if change the range

Arguments:
 - Path for the output folder

'''

from math import ceil
from math import floor
from random import randint
from random import seed
from random import shuffle
import sys
import os

def genInstance(numItems, binCapacity, l_min, l_max, instCount = 1):
  for count in range(1,instCount+1):

    lb = max(floor(binCapacity * l_min), 1)
    ub = min(ceil(binCapacity * l_max), binCapacity - 1)
    
    num = ceil(numItems / 9)
    #num = numItems
    #while ceil(num / 3) > 3:
    #  num /= 3
    
    ss = set()
    bins = []
    for i in range(int(ceil(num))):
      for x in range(100):
        i1 = randint(lb, ub)
        i2 = randint(lb, binCapacity - i1 - lb)
        i3 = binCapacity - i1 - i2
        if (i1 in ss or i2 in ss or i3 in ss) and x != 99:
          continue
        bins.append([i1, i2, i3])
        ss.add(i1)
        ss.add(i2)
        ss.add(i3)
        break
    
    while 3 * len(bins) < numItems:
      bins2 = bins.copy()
      for bin in bins2:
        bins.remove(bin)
        for i1 in bin:
          for x in range(100):
            i2 = randint(lb, binCapacity - i1 - lb)
            i3 = binCapacity - i1 - i2
            if (i2 in ss or i3 in ss) and x != 99:
              continue
            bins.append([i1, i2, i3])
            ss.add(i2)
            ss.add(i3)
            break
        
    items = []
    for bin in bins:
      items += bin
    items.sort()

    name = "N%dW%d-%03d" % (len(items), binCapacity, count)
    with open(outdirBPP + "/" + name, "w") as instFileBPP:
      instFileBPP.write("%d\n%d\n" % (len(items), binCapacity))
      for size in items:
        instFileBPP.write("%d\n" % (size))
    
  
    mp = []
    i = 0
    while i < len(items):
      mp.append([items[i], 0])
      while(i < len(items) and mp[-1][0] == items[i]):
        mp[-1][1] += 1
        i += 1

    with open(outdirCSP + "/" + name, "w") as instFileCSP:
      instFileCSP.write("%d\n%d\n" % (len(mp), binCapacity))
      for size, freq in mp:
        instFileCSP.write("%d\t%d\n" % (size, freq))
    
if __name__ == "__main__":
  if(len(sys.argv) != 2):
    print("Error: Amount of arguments is not equal to 1")
    sys.exit(1)
    
  programDir = "."

  pairs = [[200, 10000], [400, 15000], [600, 20000]]
  l_min = 0.2
  l_max = 0.4

  for numItems, binCapacity in pairs:
    infix = "_20to40-" + str(numItems)
    outdirBPP = sys.argv[1] + infix + "-BPP"
    if(not os.path.exists(outdirBPP)):
      os.mkdir(outdirBPP)
    else:
      os.system("rm " + outdirBPP + "/*")
      
    outdirCSP = sys.argv[1] + infix + "-CSP"
    if(not os.path.exists(outdirCSP)):
      os.mkdir(outdirCSP)
    else:
      os.system("rm " + outdirCSP + "/*")

    nameOutDir = sys.argv[1][sys.argv[1].rfind("/") + 1 : ]

    genInstance(numItems, binCapacity, l_min, l_max, instCount = 50)
  

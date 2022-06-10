'''
Descripttion: 
version: 
Author: zpliu
Date: 2021-02-19 21:18:32
LastEditors: zpliu
LastEditTime: 2021-02-23 22:56:39
@param: 
'''
import sys
import re


def getAllCpGsite(CpGFile):
    out = {}
    with open(CpGFile, 'r') as File:
        for line in File:
            line = line.strip("\n").split("\t")
            site = int(line[2].split("-")[-1])
            out[line[5]] = out.get(line[5], [])
            out[line[5]].append((site, line[-1]))
    return out


def readASLocation(ASevent):
    numStr = "\t".join(ASevent.split(":")[2:-1])
    # numStr = "\t".join(ASevent.split(":")[1:])
    numList = re.findall(r'\d+', numStr)
    numList = sorted([int(i) for i in numList])
    return (numList[0], numList[-1])


if __name__ == "__main__":
    ASeventLocationFile = sys.argv[1]
    CpGMethylationFile = sys.argv[2]
    out = []
    AllCpGMethlated = getAllCpGsite(CpGMethylationFile)
    with open(ASeventLocationFile, 'r') as File:
        for line in File:
            line = line.strip("\n").split("\t")
            start, end = readASLocation(line[0])
            geneId = line[0].split(";")[0]
            # print(start, end)
            # geneId = line[2]
            totalsite = 0
            DmCGsite = 0
            try:
                for site, pvalue in AllCpGMethlated[geneId]:
                    if site >= start and site <= end:
                        totalsite += 1
                        if pvalue != "None" and float(pvalue) <= 0.01:
                            DmCGsite += 1
                out.append(line[0]+"\t"+str(DmCGsite)+"\t"+str(totalsite)+"\n")
            except KeyError:
                print("nonConserved Cytosine\t"+geneId)
    with open(sys.argv[3], 'w') as File:
        for item in out:
            File.write(item)

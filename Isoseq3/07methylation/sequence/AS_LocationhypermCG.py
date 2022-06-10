'''
Descripttion: 
version: 
Author: zpliu
Date: 2021-02-19 21:18:32
LastEditors: zpliu
LastEditTime: 2021-02-22 15:45:06
@param: 
'''
import sys
import re


def getAllCpGsite(CpGFile):
    out = {}

    def getMathylation(x):
        if x == "None":
            return 0.0
        else:
            return float(x)
    with open(CpGFile, 'r') as File:
        for line in File:
            line = line.strip("\n").split("\t")
            site1 = int(line[0].split("-")[-1])
            site2 = int(line[2].split("-")[-1])
            out[line[4]] = out.get(line[4], [])
            if (getMathylation(line[6])+getMathylation(line[7])) >= 1.2:
              # hyper CG methylation site
                out[line[4]].append((site1, 1))
            if (getMathylation(line[8])+getMathylation(line[9])) >= 1.2:
              # hyper CG methylation site
                out[line[5]] = out.get(line[5], [])
                out[line[5]].append((site2, 1))
    return out


def readASLocation(ASevent):
    # numStr = "\t".join(ASevent.split(":")[2:-1])
    numStr = "\t".join(ASevent.split(":")[2:])
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
            start1, end1 = readASLocation(line[0])
            start2, end2 = readASLocation(line[3])
            # geneId = line[0].split(";")[0]
            # print(start, end)
            geneId1 = line[0].split(";")[0]
            geneId2 = line[3].split(";")[0]
            total1 = 0
            total2 = 0
            # print(start1, end1)
            try:
                for site, pvalue in AllCpGMethlated[geneId1]:
                    # non hyper mCG site
                    if site >= start1 and site <= end1:
                        total1 += 1
            except KeyError:
                print("non hyper mCG Cytosine\t"+geneId1)
            try:
                for site, pvalue in AllCpGMethlated[geneId2]:
                    if site > start2 and site < end2:
                        total2 += 1
            except KeyError:
                print("non hyper mCG Cytosine\t"+geneId2)
            out.append("\t".join(line)+"\t"+str(total1)+"\t"+str(total2)+"\n")
    with open(sys.argv[3], 'w') as File:
        for item in out:
            File.write(item)

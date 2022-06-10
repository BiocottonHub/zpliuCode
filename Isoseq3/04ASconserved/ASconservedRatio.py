'''
Descripttion: 
version: 
Author: zpliu
Date: 2021-01-06 11:00:43
LastEditors: zpliu
LastEditTime: 2021-01-06 15:27:24
@param: 
'''

import sys


def readAllASCount(ASFile):
    out = {}
    with open(ASFile, 'r') as File:
        for line in File:
            line = line.strip("\n").split("\t")
            out[line[1]] = out.get(line[1], 0)+1
    return out


def getconservedASCount(conservedASFile):
    out = {}
    with open(conservedASFile, 'r') as File:
        for line in File:
            line = line.strip("\n").split("\t")
            geneid1 = line[0].split(";")[0]
            geneid2 = line[1].split(";")[0]
            out[geneid1] = out.get(geneid1, [])
            out[geneid2] = out.get(geneid2, [])
            out[geneid1].append(line[0])
            out[geneid2].append(line[1])
    return out


if __name__ == "__main__":
    out = []
    allAScount = readAllASCount(sys.argv[1])
    conserveAScount = getconservedASCount(sys.argv[2])
    with open(sys.argv[3], 'r') as File:
        for line in File:
            line = line.strip("\n").split("\t")
            geneID1 = line[2]
            geneID2 = line[3]
            def f(x): return len(list(set(x))) if x else 0
            try:
                if allAScount[geneID1] != 0 and allAScount[geneID2] != 0:
                    conserveCount = f(conserveAScount.get(
                        geneID1, ''))+f(conserveAScount.get(geneID2, 0))
                    out.append(geneID1+"\t"+geneID2+"\t"+str(conserveCount /
                                                             (allAScount[geneID1]+allAScount[geneID2]))+"\n")
            except KeyError:
                pass
    with open(sys.argv[4], 'w') as File:
        for item in out:
            File.write(item)

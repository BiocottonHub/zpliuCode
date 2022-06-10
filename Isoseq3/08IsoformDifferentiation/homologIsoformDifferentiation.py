'''
Descripttion: 
version: 
Author: zpliu
Date: 2021-01-30 14:58:38
LastEditors: zpliu
LastEditTime: 2021-01-30 16:06:59
@param: 
'''
import sys


def getASevents(ASfile):
    out = {}
    with open(ASfile, 'r') as File:
        for line in File:
            line = line.strip("\n").split("\t")
            out[line[1]] = out.get(line[1], [])
            out[line[1]].append(
                (
                    line[2],
                    line[3].split(","),
                    line[4].split(",")
                )
            )
    return out


def getconservedIsofomList(conservedIsoformFile):
    out = {}
    with open(conservedIsoformFile, 'r') as File:
        for line in File:
            line = line.strip("\n").split("\t")
            out[line[0]] = out.get(line[0], [])
            out[line[0]].append(line[2].split("^")[1])
    return out


def getgeneIdList(geneIdFile):
    out = []
    with open(geneIdFile, 'r') as File:
        for line in File:
            line = line.strip("\n").split("\t")
            out.append(line[1])
    return out


if __name__ == "__main__":
        # all AS file
    ASeventDict = getASevents(sys.argv[1])
    conservedGroup1 = getconservedIsofomList(sys.argv[2])
    conservedGroup2 = getconservedIsofomList(sys.argv[3])
    geneIdList = getgeneIdList(sys.argv[5])
    out = []
    for geneId in geneIdList:
        for item1 in conservedGroup1[geneId]:
            for item2 in conservedGroup2[geneId]:
                try:
                    for ASevents in ASeventDict[geneId]:
                        if item1 in ASevents[1] and item2 not in ASevents[1] and item2 in ASevents[2]:
                            out.append(ASevents[0]+"\n")
                        elif item1 in ASevents[2] and item2 not in ASevents[2] and item2 in ASevents[1]:
                            out.append(ASevents[0]+"\n")
                        else:
                            pass
                except KeyError:
                    pass
    with open(sys.argv[4], 'w') as File:
        for item in list(set(out)):
            File.write(item)

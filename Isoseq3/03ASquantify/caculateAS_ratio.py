'''
Descripttion: Quantiy AS event be PacBio FL-length
version: 
Author: zpliu
Date: 2020-12-28 21:33:27
LastEditors: zpliu
LastEditTime: 2020-12-28 22:19:24
@param: 
'''
import sys
import re
#import numpy as np


def getASevent(ASeventFile: str):
    '''
    @Descripttion: 
    @param: 
    @return:  gene AS event quantify @dict
    '''
    out = {}
    with open(ASeventFile, 'r') as File:
        for line in File:
            line = line.strip("\n").split("\t")
            def f(x): return [i for i in x.split(",") if re.match(r'^PB', i)]
            supportIsoform = f(line[3])

            unsupportIsoform = [i for i in f(
                line[4]) if i not in supportIsoform]
            out[line[2]] = [supportIsoform, unsupportIsoform]
    return out


def getGeneFLread(gene_FLreadFile: str):
    '''
    @Descripttion: 
    @param: 
    @return: 
    '''
    out = {}
    with open(gene_FLreadFile, 'r') as File:
        for line in File:
            line = line.strip("\n").split("\t")
            out[line[0]] = out.get(line[0], {'count': 0})
            out[line[0]][line[1]] = int(line[5])
            out[line[0]]['count'] += int(line[5])
    return out


if __name__ == "__main__":
    ASeventDict = getASevent(sys.argv[1])
    genereadCount = getGeneFLread(sys.argv[2])
    with open(sys.argv[3], 'w') as File:
        for event, IsoformArray in ASeventDict.items():
            if len(IsoformArray[1]) == 0:
                # PacBio compared with reference transcript so skip
                pass
            else:
                geneId = event.split(";")[0]
                supportreadCount = 0
                unsupportreadCount = 0
                for i in IsoformArray[0]:
                    supportreadCount += genereadCount[geneId][i]
                for i in IsoformArray[1]:
                    unsupportreadCount += genereadCount[geneId][i]
                File.write(geneId+"\t"+event+"\t"+str(supportreadCount) +
                           "\t"+str(unsupportreadCount)+"\t"+str(genereadCount[geneId]['count'])+"\n")

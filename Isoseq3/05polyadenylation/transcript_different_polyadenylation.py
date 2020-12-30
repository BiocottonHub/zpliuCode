'''
Descripttion: enumerate the different isoform to comapare
version: 1.0
Author: zpliu
Date: 2020-12-30 15:43:15
LastEditors: zpliu
LastEditTime: 2020-12-30 19:24:00
@param: 
'''
from itertools import combinations
from collections import Counter
import sys
sys.path.insert(0, '/public/home/zpliu/github/zpliuCode/script/genestruct')

while 1:
    #from readgtf import getTranscriptInfo
    from readgtf import getGeneInfo
    break


def deletTranscript(isoform1Object, isoform2Object):
    exon1 = sorted(isoform1Object.getExonCoordinate())
    exon2 = sorted(isoform2Object.getExonCoordinate())
    isoformName1 = isoform1Object.transcriptName
    isoformName2 = isoform2Object.transcriptName
    stand = isoform1Object.stand
    geneId = isoform2Object.geneId
    countTime = Counter(exon1+exon2)
    plotASame = ''
    TSSSame = ''
    splitsitSame = ''
    # All splice site same
    tmp = [True if countTime[i] == 2 else False for i in exon1[1:-1]]
    # TSS site
    if stand == "+" and countTime[exon1[0]] == 2:
        TSSSame = 'TSSsame'
    if stand == "+" and countTime[exon1[0]] == 1:
        TSSSame = 'TSSdiff'
    if stand == "-" and countTime[exon1[-1]] == 2:
        TSSSame = 'TSSsame'
    if stand == "-" and countTime[exon1[-1]] == 1:
        TSSSame = 'TSSdiff'
    # ployA site
    if stand == "-" and countTime[exon1[0]] == 2:
        plotASame = 'ployAsame'
    if stand == "-" and countTime[exon1[0]] == 1:
        plotASame = 'ployAdiff'
    if stand == "+" and countTime[exon1[-1]] == 2:
        plotASame = 'ployAsame'
    if stand == "+" and countTime[exon1[-1]] == 1:
        plotASame = 'ployAdiff'
    if len(exon1) != len(exon2):
        splitsitSame = 'splitDiff'
    elif all(tmp):
        splitsitSame = 'splitsame'
    else:
        splitsitSame = 'splitDiff'

    return "\t".join([geneId, isoformName1, isoformName2, TSSSame, plotASame, splitsitSame])+"\n"


if __name__ == "__main__":
    gtfFile = sys.argv[1]
    outFile = sys.argv[2]
    geneMessage = getGeneInfo(gtfFile)
    with open(outFile, 'w') as File:
        for key in geneMessage:
            if len(geneMessage[key].transcripts) > 1:
                for isoform1, isoform2 in combinations([i for i in geneMessage[key].transcriptsObject], 2):
                    result = deletTranscript(isoform1, isoform2)
                    # print(result)
                    File.write(result)
            else:
                pass

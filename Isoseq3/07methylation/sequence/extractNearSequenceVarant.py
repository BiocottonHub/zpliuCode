'''
Descripttion: 
version: 
Author: zpliu
Date: 2021-02-21 16:34:09
LastEditors: zpliu
LastEditTime: 2021-02-21 17:20:32
@param: 
'''
import sys
from extractASsequence import getsequenceByLoction


def getVarianLocation(evetLocation, varianSite, strand):
    # deal with donor file
    chromsome = evetLocation.split(":")[0]
    start = int(evetLocation.split(":")[1].split("-")[0])
    end = int(evetLocation.split("-")[1])
    if strand == "+":
        site = start+varianSite-1
        return chromsome+":"+str(site-10)+"-"+str(site+10)
    else:
        site = end-varianSite+1
        return chromsome+":"+str(site-10)+"-"+str(site+10)


def readVariantSite(variantFile, genomeSequenceFile1, genomeSequenceFile2):
    out = []
    with open(variantFile, 'r') as File:
        for line in File:
            line = line.strip("\n").split("\t")
            site1 = line[4]
            site2 = line[5]
            strand1 = line[1]
            strand2 = line[3]
            Location1 = getVarianLocation(site1, int(line[6]), strand1)
            Location2 = getVarianLocation(site2, int(line[6]), strand2)
            sequence1 = getsequenceByLoction(
                genomeSequenceFile1, Location1, strand1)
            sequence2 = getsequenceByLoction(
                genomeSequenceFile2, Location2, strand2)
            out.append(sequence1+"\n")
            out.append(sequence2+"\n")
    return out


if __name__ == "__main__":
    genomeFile1 = sys.argv[1]
    genomeFile2 = sys.argv[2]
    varianSiteFile = sys.argv[3]
    outFile = sys.argv[4]
    with open(outFile, 'w') as File:
        for item in readVariantSite(varianSiteFile, genomeFile1, genomeFile2):
            File.write(item)

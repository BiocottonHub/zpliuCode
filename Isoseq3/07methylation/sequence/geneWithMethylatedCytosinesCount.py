'''
Descripttion: 
version: 
Author: zpliu
Date: 2021-02-18 19:48:13
LastEditors: zpliu
LastEditTime: 2021-02-19 11:04:34
@param: 
'''
import sys
import os
from extractASsequence import getsequenceByLoction


def getMethlatedCytosineCount(MethylationFile, Location):

    MethlationLines = os.popen(
        'samtools view -O SAM %s %s' % (MethylationFile, Location)).read()
    if MethlationLines == "":
        return '0'
    else:
        count = 0
        for item in MethlationLines.split("\n")[:-1]:
            read1 = item.split("\t")[0].split("-")[2]
            read2 = item.split("\t")[0].split("-")[3]
            if eval(read1+"+"+read2) >= 3:
                count += 1
        return str(count)


def getCGCount(sequence, strand):
    return str(sequence.count("CG"))


if __name__ == "__main__":

    MethylationFile1 = sys.argv[1]
    MethylationFile2 = sys.argv[2]
    LocationFile = sys.argv[3]
    genomeSequence = sys.argv[4]
    out = []
    with open(LocationFile, 'r') as File:
        File.readline()
        for line in File:
            line = line.strip("\n").split("\t")
            location = line[0]
            strand = line[2]
            sequence = getsequenceByLoction(genomeSequence, location, strand)
            methlatedCount1 = getMethlatedCytosineCount(
                MethylationFile1, location)
            methlatedCount2 = getMethlatedCytosineCount(
                MethylationFile2, location)
            totalCGCount = getCGCount(sequence, strand)
            out.append(location+"\t"+line[1]+"\t"+line[2]+"\t" + totalCGCount+"\t" +
                       methlatedCount1+"\t"+methlatedCount2+"\n")
    with open(sys.argv[3], 'w') as File:
        File.write("#location\t" +
                   "geneId\tstrand\ttotalCGCount\trep1\trep2\n")
        for item in out:
            File.write(item)

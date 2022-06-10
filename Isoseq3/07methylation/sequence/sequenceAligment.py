'''
Descripttion:
version:
Author: zpliu
Date: 2021-02-14 16:54:56
LastEditors: zpliu
LastEditTime: 2021-02-21 16:27:31
@param:
'''
from extractASsequence import getsequenceByLoction
import sys
import os
import re
from itertools import groupby


def getLocations(locationFile):
    out = []
    with open(locationFile, 'r') as File:
        for line in File:
            line = line.strip("\n").split("\t")
            out.append(line)
    return out


def getDonorSequenceVariant(muscleOutStr, location1, location2, geneId1, geneId2, strand1, strand2):
    out = []
    index = 0
    tmp = []
    for flagment in muscleOutStr.split("\n\n")[1::]:
        sequenceflagment1, sequenceflagment2 = [
            re.split("\s+", i)[-1] for i in flagment.strip("\n").split("\n")[0:2]]
        conserveflagment = flagment.strip("\n").split(
            "\n")[2][::-1][:len(sequenceflagment1)][::-1]
        print(conserveflagment)
        for sequence1, sequence2, conserved in zip(sequenceflagment1, sequenceflagment2, conserveflagment):
            if conserved == "*":
                index += 1
            else:
                if re.search(r'([ATCGN][ATCGN])', sequence1+sequence2):
                    index += 1
                    out.append(geneId1+"\t"+strand1+"\t"+geneId2+"\t"+strand2+"\t"+location1+"\t"+location2 +
                               "\t"+str(index)+"\t1\tSNP")
                else:
                    index += 1
                    tmp.append(index)
    # 获取连续的缺失范围及分组
    # print(out)
    # print(tmp)
    def fun(x): return x[1]-x[0]
    for k, g in groupby(enumerate(tmp), fun):
        l1 = [j for i, j in g]
        # print(l1)
        if len(l1) <= 50:
            out.append(geneId1+"\t"+strand1+"\t"+geneId2+"\t"+strand2+"\t"+location1+"\t"+location2 +
                       "\t"+str(l1[0])+"\t"+str(len(l1))+"\tindel")
        else:
            out.append(geneId1+"\t"+strand1+"\t"+geneId2+"\t"+strand2+"\t"+location1+"\t"+location2 +
                       "\t"+str(l1[0])+"\t"+str(len(l1))+"\tSV")

        # for item in l1:
        #     print(l1)

    return out


def getAcceptorSequenceVariant(muscleOutStr, location1, location2, geneId1, geneId2, strand1, strand2):
    out = []
    index = 0
    tmp = []
    for flagment in muscleOutStr.split("\n\n")[::-1][0:-1]:
        sequenceflagment1, sequenceflagment2 = [
            re.split("\s+", i)[-1][::-1] for i in flagment.strip("\n").split("\n")[0:2]]
        conserveflagment = flagment.strip("\n").split(
            "\n")[2][::-1][:len(sequenceflagment1)]
        # print(conserveflagment)
        for sequence1, sequence2, conserved in zip(sequenceflagment1, sequenceflagment2, conserveflagment):
            if conserved == "*":
                index += 1
            else:
                if re.search(r'([ATCGN][ATCGN])', sequence1+sequence2):
                    index += 1
                    out.append(geneId1+"\t"+strand1+"\t"+geneId2+"\t"+strand2+"\t"+location1+"\t"+location2 +
                               "\t"+str(index)+"\t1\tSNP")
                else:
                    index += 1
                    tmp.append(index)
    # 获取连续的缺失范围及分组
    # print(out)
    # print(tmp)
    def fun(x): return x[1]-x[0]
    for k, g in groupby(enumerate(tmp), fun):
        l1 = [j for i, j in g]
        if len(l1) <= 50:
            out.append(geneId1+"\t"+strand1+"\t"+geneId2+"\t"+strand2+"\t"+location1+"\t"+location2 +
                       "\t"+str(l1[0])+"\t"+str(len(l1))+"\tindel")
        else:
            out.append(geneId1+"\t"+strand1+"\t"+geneId2+"\t"+strand2+"\t"+location1+"\t"+location2 +
                       "\t"+str(l1[0])+"\t"+str(len(l1))+"\tSV")

        # for item in l1:
        #     print(l1)

    return out


def sequenceAligment(genomeSequenceFile1, genomeSequenceFile2, locations1, locations2, DonroOutFile, AcceptorOutFile):
    out = []
    for line1, line2 in zip(locations1, locations2):
        location1 = line1[0]+":" + \
            str(int(line1[1])-51)+"-" + \
            str(int(line1[2])+51)  # contained 20bp exon region
        location2 = line2[0]+":" + \
            str(int(line2[1])-51)+"-" + \
            str(int(line2[2])+51)  # contained 20bp exon region
        strand1 = line1[4]
        strand2 = line2[4]
        geneId1 = line1[3].split(";")[0]
        geneId2 = line2[3].split(";")[0]
        sequence1 = getsequenceByLoction(
            genomeSequenceFile1, location1, strand1)
        sequence2 = getsequenceByLoction(
            genomeSequenceFile2, location2, strand2)
        infasta = sequence1+"\n"+sequence2
        tmp = os.popen(
            "printf \"%s\"|~/software/muscle3.8.31_i86linux64 -clw  2>/dev/null" % infasta).read()
        eventLength1 = eval(line1[2]+"-"+line1[1])+1
        eventLength2 = eval(line2[2]+"-"+line2[1])+1
        print(tmp)
        # print(getSequenceVariant(tmp, location1, location2))
        with open(DonroOutFile, 'a') as File:
            for item in getDonorSequenceVariant(tmp, location1, location2, geneId1, geneId2, strand1, strand2):
                File.write(item+"\n")
        with open(AcceptorOutFile, 'a') as File:
            for item in getAcceptorSequenceVariant(tmp, location1, location2, geneId1, geneId2, strand1, strand2):
                File.write(item+"\n")
        out.append(location1+"\t"+strand1+"\t" + geneId1+"\t" +
                   location2+"\t"+strand2 + "\t" + geneId2+"\t"+str(eventLength1)+"\t" +
                   str(eventLength2)+"\t" +
                   str(len(re.findall(r'\*', tmp))))
    return out


if __name__ == "__main__":
    genomeSequenceFile1 = sys.argv[1]
    genomeSequenceFile2 = sys.argv[2]
    eventLocaltion1 = getLocations(sys.argv[3])
    eventLocaltion2 = getLocations(sys.argv[4])
    DonroOutFile = sys.argv[5]
    AcceptorOutFile = sys.argv[6]
    with open(sys.argv[5], 'w') as File:
        pass
    with open(sys.argv[6], 'w') as File:
        pass
    out = sequenceAligment(genomeSequenceFile1, genomeSequenceFile2,
                           eventLocaltion1, eventLocaltion2, DonroOutFile, AcceptorOutFile)
    with open(sys.argv[7], 'w') as File:
        for item in out:
            File.write(item+"\n")

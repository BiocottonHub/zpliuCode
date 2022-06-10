'''
Descripttion: 
version: 
Author: zpliu
Date: 2021-02-18 15:59:58
LastEditors: zpliu
LastEditTime: 2021-02-18 21:23:59
@param: 
'''
from extractASsequence import getsequenceByLoction
import sys
import os
import re


def getLocations(locationFile):
    out = []
    with open(locationFile, 'r') as File:
        for line in File:
            line = line.strip("\n").split("\t")
            out.append(line)
    return out


def geconservedCytosine(muscleOutStr, location1, location2, strand1, strand2):
    out = []
    index1 = 0
    index2 = 0
    for flagment in muscleOutStr.split("\n\n")[1::]:
        sequenceflagment1, sequenceflagment2 = [
            re.split("\s+", i)[-1] for i in flagment.strip("\n").split("\n")[0:2]]
        conserveflagment = flagment.strip("\n").split(
            "\n")[2][::-1][:len(sequenceflagment1)][::-1]
        for sequence1, sequence2, conserved in zip(sequenceflagment1, sequenceflagment2, conserveflagment):
            if conserved == "*" and re.search(r'([Cc][Cc])', sequence1+sequence2):
                index1 += 1
                index2 += 1
                # print(sequence1+sequence2+conserved)
                out.append((index1, index2))
            elif re.search(r'([-][ATCGatcgNn])', sequence1+sequence2):
                index2 += 1
            elif re.search(r'([ATCGatcgNn][-])', sequence1+sequence2):
                index1 += 1
            else:
                index1 += 1
                index2 += 1
    # print(out)
    out2 = []
    if strand1 == "+" and strand2 == "+":
        start1 = int(location1[1])
        start2 = int(location2[1])
        for index1, index2 in out:
            out2.append(
                location1[0]+":"+str(start1+index1-1)+"-"+str(start1+index1-1)+"\t"+strand1+"\t" +
                location2[0]+":"+str(start2+index2-1)+"-" +
                str(start2+index2-1)+"\t"+strand2
            )
    elif strand1 == "-" and strand2 == "-":
        start1 = int(location1[2])
        start2 = int(location2[2])
        for index1, index2 in out:
            out2.append(
                location1[0]+":"+str(start1-index1+1)+"-"+str(start1-index1+1)+"\t"+strand1+"\t" +
                location2[0]+":"+str(start2-index2+1)+"-"+str(start2-index2+1)
                + "\t"+strand2
            )

    elif strand1 == "-" and strand2 == "+":
        start1 = int(location1[2])
        start2 = int(location2[1])
        for index1, index2 in out:
            out2.append(
                location1[0]+":"+str(start1-index1+1)+"-"+str(start1-index1+1)+"\t"+strand1+"\t" +
                location2[0]+":"+str(start2+index2-1)+"-" +
                str(start2+index2-1)+"\t"+strand2
            )

    elif strand1 == "+" and strand2 == "-":
        start1 = int(location1[1])
        start2 = int(location2[2])
        for index1, index2 in out:
            out2.append(
                location1[0]+":"+str(start1+index1-1)+"-"+str(start1+index1-1)+"\t"+strand1+"\t" +
                location2[0]+":"+str(start2-index2+1)+"-"+str(start2-index2+1)
                + "\t"+strand2
            )
    return out2


def sequenceAligment(genomeSequenceFile1, genomeSequenceFile2, locations1: list, locations2, consevedCytosineFile):
    out = []
    for line1, line2 in zip(locations1, locations2):
        location1 = line1[0]+":" + \
            str(int(line1[1]))+"-" + \
            str(int(line1[2]))  # contained 20bp exon region
        location2 = line2[0]+":" + \
            str(int(line2[1]))+"-" + \
            str(int(line2[2]))  # contained 20bp exon region
        strand1 = line1[5]
        strand2 = line2[5]
        sequence1 = getsequenceByLoction(
            genomeSequenceFile1, location1, strand1)
        sequence2 = getsequenceByLoction(
            genomeSequenceFile2, location2, strand2)
        infasta = sequence1+"\n"+sequence2
        # sequence too long so except OSError: [Errno 7] Argument list too long
        try:
            tmp = os.popen(
                "printf \"%s\"|~/software/muscle3.8.31_i86linux64 -clw  2>/dev/null" % infasta).read()
            out = geconservedCytosine(tmp, line1, line2, strand1, strand2)
            with open(consevedCytosineFile, 'a') as File:
                for item in out:
                    File.write(item+"\t"+line1[3]+"\t"+line2[3]+"\n")
        except OSError as error:
            # sequence too long
            print(error)
            print(location1)
            pass
        except UnicodeDecodeError as error:
            # muscle error
            print(error)
            print(location1)
            pass


if __name__ == "__main__":
    genomeFile1 = sys.argv[1]
    genomeFile2 = sys.argv[2]
    locationList1 = getLocations(sys.argv[3])
    locationList2 = getLocations(sys.argv[4])
    with open(sys.argv[5], 'w') as File:
        pass
    sequenceAligment(genomeFile1, genomeFile2, locationList1,
                     locationList2, sys.argv[5])

'''
Descripttion: 
version: 
Author: zpliu
Date: 2021-02-07 13:20:12
LastEditors: zpliu
LastEditTime: 2021-02-09 17:18:25
@param: 
'''
import sys
import os
import re


def getCGmethylation(OutDict, genomeSequence, MethylationFile, eventCoordinate):
    flagmentLength = eventCoordinate[2]-eventCoordinate[1]
    if eventCoordinate[-1] == "+" and flagmentLength >= 74:
        for index, start, end in zip(range(1, 276, 1), range(eventCoordinate[1]-200, eventCoordinate[1]+75), range(eventCoordinate[1]-199, eventCoordinate[1]+76)):
            coordinate = eventCoordinate[0]+":"+str(start)+"-"+str(end)
            delwithDinucle(OutDict, index, coordinate,
                           genomeSequence, MethylationFile, eventCoordinate[-1])
        for index, start, end in zip(range(276, 551, 1), range(eventCoordinate[2]-75, eventCoordinate[2]+200), range(eventCoordinate[2]-74, eventCoordinate[2]+201)):
            coordinate = eventCoordinate[0]+":"+str(start)+"-"+str(end)
            delwithDinucle(OutDict, index, coordinate,
                           genomeSequence, MethylationFile, eventCoordinate[-1])
    if eventCoordinate[-1] == "+" and flagmentLength < 74:
        for index, start, end in zip(range(1, 200+flagmentLength+2, 1), range(eventCoordinate[1]-200, eventCoordinate[1]+flagmentLength+1), range(eventCoordinate[1]-199, eventCoordinate[1]+flagmentLength+1)):
            coordinate = eventCoordinate[0]+":"+str(start)+"-"+str(end)
            delwithDinucle(OutDict, index, coordinate,
                           genomeSequence, MethylationFile, eventCoordinate[-1])
        for index, start, end in zip(range(350-flagmentLength-1, 551, 1), range(eventCoordinate[2]-flagmentLength-1, eventCoordinate[2]+200), range(eventCoordinate[2]-flagmentLength, eventCoordinate[2]+201)):
            coordinate = eventCoordinate[0]+":"+str(start)+"-"+str(end)
            delwithDinucle(OutDict, index, coordinate,
                           genomeSequence, MethylationFile, eventCoordinate[-1])
    if eventCoordinate[-1] == "-" and flagmentLength >= 74:
        for index, start, end in zip(range(1, 276, 1), range(eventCoordinate[2]+200, eventCoordinate[2]-75, -1), range(eventCoordinate[2]+201, eventCoordinate[2]-74, -1)):
            coordinate = eventCoordinate[0]+":"+str(start)+"-"+str(end)
            delwithDinucle(OutDict, index, coordinate,
                           genomeSequence, MethylationFile, eventCoordinate[-1])
        for index, start, end in zip(range(276, 551, 1), range(eventCoordinate[1]+75, eventCoordinate[1]-200, -1), range(eventCoordinate[1]+76, eventCoordinate[1]-199, -1)):
            coordinate = eventCoordinate[0]+":"+str(start)+"-"+str(end)
            delwithDinucle(OutDict, index, coordinate,
                           genomeSequence, MethylationFile, eventCoordinate[-1])
    if eventCoordinate[-1] == "-" and flagmentLength < 74:
        for index, start, end in zip(range(1, 200+flagmentLength+2, 1), range(eventCoordinate[2]+200, eventCoordinate[2]-flagmentLength-1, -1), range(eventCoordinate[2]+201, eventCoordinate[2]-flagmentLength, -1)):
            coordinate = eventCoordinate[0]+":"+str(start)+"-"+str(end)
            delwithDinucle(OutDict, index, coordinate,
                           genomeSequence, MethylationFile, eventCoordinate[-1])
        for index, start, end in zip(range(350-flagmentLength-1, 551, 1), range(eventCoordinate[1]+flagmentLength+1, eventCoordinate[1]-200, -1), range(eventCoordinate[1]+flagmentLength+2, eventCoordinate[1]-199, -1)):
            coordinate = eventCoordinate[0]+":"+str(start)+"-"+str(end)
            delwithDinucle(OutDict, index, coordinate,
                           genomeSequence, MethylationFile, eventCoordinate[-1])


def delwithDinucle(outDict, index, coordinate, genomeSequence, MethylationFile, stand):
    sequence = os.popen("samtools faidx %s %s" %
                        (genomeSequence, coordinate)).read()
    sequence = sequence.split("\n")[1]
    if re.match('CG', sequence, re.IGNORECASE) and stand == "+":
        methylationline = os.popen(
            "samtools view -O SAM  %s %s" % (MethylationFile, coordinate)).read()
        if methylationline != '' and float(methylationline.split("\t")[0].split("-")[-1]) <= 0.1:
            # GC count but did not detect methylation
            outDict[index][0] += 1
            outDict[index][1] += 1
        elif methylationline != '' and float(methylationline.split("\t")[0].split("-")[-1]) > 0.1:
            outDict[index][1] += 1
    elif re.match('GC', sequence, re.IGNORECASE) and stand == "-":
        methylationline = os.popen(
            "samtools view -O SAM  %s %s" % (MethylationFile, coordinate)).read()
        if methylationline != '' and float(methylationline.split("\t")[0].split("-")[-1]) <= 0.1:
            # GC count but did not detect methylation
            outDict[index][0] += 1
            outDict[index][1] += 1
        elif methylationline != '' and float(methylationline.split("\t")[0].split("-")[-1]) > 0.1:
            outDict[index][1] += 1
    else:
        pass


def initializationDict():
    out = {}
    for i in range(1, 551):
        out[i] = out.get(i, [0, 0])  # methylation count ,CG count
    return out


if __name__ == "__main__":
    OutDict = initializationDict()
    genomeSequenceFile = sys.argv[1]
    CpGmethylationFile = sys.argv[2]
    eventCooridinateFile = sys.argv[3]
    outFile = sys.argv[4]
    with open(eventCooridinateFile, 'r') as File:
        for line in File:
            line = line.strip("\n").split("\t")
            line[1] = int(line[1])
            line[2] = int(line[2])
            getCGmethylation(OutDict, genomeSequenceFile,
                             CpGmethylationFile, line)
    with open(outFile, 'w') as File:
        for key, values in OutDict.items():
            File.write(str(key)+"\t"+str(values[0])+"\t"+str(values[1])+"\n")

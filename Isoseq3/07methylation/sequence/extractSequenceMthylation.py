'''
Descripttion:
version:
Author: zpliu
Date: 2021-02-16 11:25:06
LastEditors: zpliu
LastEditTime: 2021-02-16 11:36:50
@param:
'''
import sys
import os


def getMethylationRatio(samtoolsOut):
    if samtoolsOut == "":
        return "0\t0"
    else:
        supportMethCount = 0
        unsupportMethCount = 0
        for item in samtoolsOut.split("\n")[:-1]:
            supportMethCount += int(item.split("-")[2])
            unsupportMethCount += int(item.split("-")[3])
        return str(supportMethCount)+"\t"+str(unsupportMethCount)


if __name__ == "__main__":
    genomeMethylation1 = sys.argv[1]
    genomeMethylation2 = sys.argv[2]
    locationFile = sys.argv[3]
    out = []
    with open(locationFile, 'r') as File:
        for line in File:
            line = line.strip("\n").split("\t")
            MethlationRatio1 = getMethylationRatio(os.popen(
                'samtools view -O SAM %s %s' % (genomeMethylation1, line[0])).read())
            MethlationRatio2 = getMethylationRatio(os.popen(
                'samtools view -O SAM %s %s' % (genomeMethylation2, line[2])).read())
            out.append(line[0]+"\t"+line[2]+"\t" +
                       MethlationRatio1+"\t"+MethlationRatio2+"\n")
    with open(sys.argv[4], 'w') as File:
        for item in out:
            File.write(item)

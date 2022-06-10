'''
Descripttion: 
version: 
Author: zpliu
Date: 2021-02-07 19:20:06
LastEditors: zpliu
LastEditTime: 2021-02-13 09:30:29
@param: 
'''
import sys
import os
import re


def getRegionMethylation(RegionFile, MethylationFile):
    out = {}
    with open(RegionFile, 'r') as File:
        for line in File:
            line = line.strip("\n").split("\t")
            Region = line[0]+":"+line[1]+"-"+line[2]
            MethlationLines = os.popen(
                'samtools view -O SAM %s %s' % (MethylationFile, Region)).read()
            if MethlationLines == "":
                out[line[3]] = out.get(line[3], '0\t0\t0')
            else:
                supportMethCount = 0
                unsupportMethCount = 0
                for item in MethlationLines.split("\n")[:-1]:
                    pvalue = item.split("\t")[0].split("-")[4]
                    if re.search('e', pvalue):
                        supportMethCount += 1
                    elif float(pvalue) <= 0.01:
                        supportMethCount += 1
                    else:
                        unsupportMethCount += 1
                out[line[3]] = out.get(line[3], str(
                    supportMethCount/(supportMethCount+unsupportMethCount)))
    return out


if __name__ == "__main__":
    regionFile = sys.argv[1]
    methylationFile = sys.argv[2]
    outFile = sys.argv[3]
    with open(outFile, 'w') as File:
        for key, value in getRegionMethylation(regionFile, methylationFile).items():
            File.write(key+"\t"+value+"\n")

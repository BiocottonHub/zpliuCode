'''
Descripttion: 
version: 
Author: zpliu
Date: 2020-12-26 20:26:15
LastEditors: zpliu
LastEditTime: 2021-01-24 09:36:48
@param: 
'''
import random
import sys
import numpy as np
'''
+ All gene containde full read count File
+ Full reads count
'''


def readFLreadCount(FLreadCountFile: str):
    '''
    @Descripttion: name each Fl-reads to the dict
    @param: 
    @return: 
    '''
    out = {}
    count = 1
    with open(FLreadCountFile, 'r') as File:
        for line in File:
            line = line.strip("\n").split("\t")
            for i in range(count, count+int(line[2]), 1):
                out[i] = [line[0], line[1]]
            count += int(line[2])
    return out


if __name__ == "__main__":
    totalreadCount = int(sys.argv[2])
    readSamples = readFLreadCount(sys.argv[1])
    geneCount = []
    isoformCount = []
    out = []
    for i in np.arange(0.01, 1.01, 0.01):
        sampleCount = round(totalreadCount*i)
        for readID in random.sample(range(1, totalreadCount+2, 1), sampleCount):
            try:
                geneCount.append(readSamples[readID][0])
                isoformCount.append(readSamples[readID][1])
            except KeyError:
                pass  # ignore the full read mapping to Intergenic Region
        out.append("\t".join(
            [str(i), str(len(list(set(geneCount)))), str(len(set(isoformCount)))])+"\n")
        geneCount = []
        isoformCount = []
    with open(sys.argv[3], 'w') as File:
        for item in out:
            File.write(item)

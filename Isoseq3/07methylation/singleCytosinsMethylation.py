'''
Descripttion:
version:
Author: zpliu
Date: 2021-01-08 11:14:12
LastEditors: zpliu
LastEditTime: 2021-01-11 19:14:19
@param:  some error
'''

import re
import sys


def readMethylationFile(methylationFile):
    out = {}
    with open(methylationFile, 'r') as File:
        for line in File:
            line = line.strip("\n").split("\t")
            key = line[0]+"-"+line[1]
            out[key] = line[4]+"\t"+line[5]
    return out


def readConservedCytosinsFile(conservedCytosinFile):
    with open(conservedCytosinFile, 'r') as File:
        for line in File:
            line = line.strip("\n").split("\t")
            chromsome1 = line[4].split(":")[1]
            chromsome2 = line[5].split(":")[1]
            if re.match('^-', line[0]):
                yield (None, chromsome2+"-"+line[3], line[4], line[5])
            elif re.match('^-', line[1]):
                yield (chromsome2+"-"+line[2], None, line[4], line[5])
            else:
                yield (chromsome1+"-"+line[2], chromsome2+"-"+line[3], line[4], line[5])


def writeFile(outFile, outItems):
    with open(outFile, 'a+') as File:
        for item in outItems:
            File.write(item+"\n")


if __name__ == "__main__":
    methtylaTionDict1 = readMethylationFile(sys.argv[1])
    methtylaTionDict2 = readMethylationFile(sys.argv[2])
    with open(sys.argv[4], 'w') as File:
        pass
    for key1, key2, event1, event2 in readConservedCytosinsFile(sys.argv[3]):
        try:
            item1 = methtylaTionDict1[key1]
        except KeyError:
            item1 = 'None\tNone'
        try:
            item2 = methtylaTionDict2[key2]
        except KeyError:
            item2 = 'None\tNone'
        item = (
            "\t".join([str(key1), str(key2),
                       item1, item2,
                       event1, event2]),)
        writeFile(sys.argv[4], item)

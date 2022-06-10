'''
Descripttion: 
version: 
Author: zpliu
Date: 2021-01-05 21:11:00
LastEditors: zpliu
LastEditTime: 2021-01-05 22:03:12
@param: 
'''
import re
import sys


def getmostSameAS(ASlength1, ASArray):
    global ASlength
    tmp = 999999999
    event = ''
    for AS in ASArray:
        if abs(ASlength1-ASlength[AS]) < tmp:
            tmp = abs(ASlength1-ASlength[AS])
            event = AS
    return event


if __name__ == "__main__":
    out = []
    ASDict = {}
    ASlength = {}
    with open(sys.argv[1], 'r') as File:
        for line in File:
            line = line.strip("\n").split("\t")
            ASDict[line[0]] = ASDict.get(line[0], [])
            ASDict[line[0]].append(line[1])
            ASDict[line[1]] = ASDict.get(line[1], [])
            ASDict[line[1]].append(line[0])
            ASlength[line[0]] = int(line[-1])
            ASlength[line[1]] = int(line[-2])
    for key1 in ASDict:
        for item in ASDict[key1]:
            if item == key1:
                pass
            else:
                for key2 in ASDict[item]:
                    if key2 not in ASDict[key1]:
                        ASDict[key1].append(key2)
                        for key3 in ASDict[key2]:
                            if key3 not in ASDict[key1]:
                                ASDict[key1].append(key3)
                                for key4 in ASDict[key3]:
                                    if key4 not in ASDict[key1]:
                                        ASDict[key1].append(key4)
    with open(sys.argv[2], 'w') as File:
        for key in ASDict:
            A2AS = [i for i in ASDict[key] if re.match(r'evm', i)]
            AtAS = [i for i in ASDict[key] if re.match(r'Ghir_A', i)]
            DtAS = [i for i in ASDict[key] if re.match(r'Ghir_D', i)]
            D5AS = [i for i in ASDict[key] if re.match(r'Gor', i)]
            tmp2 = []
            tmp2.append(getmostSameAS(ASlength[key], A2AS))
            tmp2.append(getmostSameAS(ASlength[key], AtAS))
            tmp2.append(getmostSameAS(ASlength[key], DtAS))
            tmp2.append(getmostSameAS(ASlength[key], D5AS))
            File.write("\t".join(tmp2)+"\n")

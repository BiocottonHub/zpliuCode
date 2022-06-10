'''
Descripttion: 
version: 
Author: zpliu
Date: 2021-01-05 21:11:00
LastEditors: zpliu
LastEditTime: 2021-01-06 18:57:31
@param: 
'''
import re
import sys

'''
##获得同源基因间并行分化与趋同进化
'''

if __name__ == "__main__":
    out = []
    ASDict = {}

    with open(sys.argv[1], 'r') as File:
        for line in File:
            line = line.strip("\n").split("\t")
            ASDict[line[0]] = ASDict.get(line[0], [])
            ASDict[line[0]].append(line[1])
            ASDict[line[1]] = ASDict.get(line[1], [])
            ASDict[line[1]].append(line[0])

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
            tmp = sorted(list(set([key]+ASDict[key])))
            File.write("\t".join(tmp)+"\n")

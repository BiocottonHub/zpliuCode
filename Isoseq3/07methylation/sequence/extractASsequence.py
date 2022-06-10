'''
Descripttion: 
version: 
Author: zpliu
Date: 2021-02-14 16:39:53
LastEditors: zpliu
LastEditTime: 2021-02-17 16:07:13
@param: 
'''
import sys
import os


def reverseSequence(sequence):
    dict1 = {
        'A': 'T',
        'T': 'A',
        'G': 'C',
        'C': 'G',
        'N': 'N',
        'a': 'T',
        't': 'A',
        'g': 'C',
        'c': 'G',
        'n': 'N'
    }
    item = ''
    for i in sequence:
        item += dict1[i]
    return item[::-1]


def getsequenceByLoction(genomeSeuenceFile, location, strand):
    sequence = os.popen("samtools faidx % s % s" %
                        (genomeSeuenceFile, location)).read()
    id = ''
    Anothersequence = ''
    for item in sequence.split("\n"):
        if item.startswith(">"):
            id = item
        else:
            Anothersequence += item
    if strand == '+':
        return id+"\n"+Anothersequence
    else:
        return id+"\n"+reverseSequence(Anothersequence)


if __name__ == "__main__":
    genomeSeuenceFile = sys.argv[1]
    locaitonFile = sys.argv[2]
    outFile = sys.argv[3]
    out = []
    with open(locaitonFile, 'r') as File:
        for line in File:
            line = line.strip("\n").split("\t")
            location = line[0]+":" + \
                str(int(line[1])-50)+"-"+str(int(line[2])+50)
            out.append(getsequenceByLoction(
                genomeSeuenceFile, location, line[4])+"\n")
    with open(outFile, 'w') as File:
        for item in out:
            File.write(item)

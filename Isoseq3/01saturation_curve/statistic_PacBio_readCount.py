'''
Descripttion: 
version: 
Author: zpliu
Date: 2020-12-21 09:48:17
LastEditors: zpliu
LastEditTime: 2020-12-26 19:48:49
@param: 
'''
import re
import argparse
import numpy as np

'''
+ quiver: total full length fasta file 
+ collapse: collapse isoform into transcripts file 
+ gene: contained gene and PacBio ID the sparate by \t
+ out: out file ,each line  with full length and non-full length read count
'''


def getROIReads(quiverFile: str):
    """
    docstring
    """
    readCount = {}
    pattern1 = r'>([^ ]*)'
    pattern2 = r'coverage=([0-9]*)'
    with open(quiverFile, 'r') as File:
        for line in File:
            if re.match(r'^>', line):
                id = re.findall(pattern1, line)[0]
                count = re.findall(pattern2, line)
                readCount[id] = np.array([int(i) for i in count])
    return readCount


def getquiverId(collapseFile):
    """
    docstring
    """
    collapse = {}
    with open(collapseFile, 'r') as File:
        for line in File:
            line = line.strip("\n").split("\t")
            collapse[line[0]] = line[1].split(",")
    return collapse


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-quiver', help='polished reads')
    parser.add_argument('-collapse', help='isoform group reads')
    parser.add_argument('-gene', help='geneId')
    parser.add_argument('-out', help='out file name')
    namespace = parser.parse_args()
    out = []
    collapse = getquiverId(namespace.collapse)
    ROIreads = getROIReads(namespace.quiver)
    with open(namespace.gene, 'r') as File:
        for line in File:
            line = line.strip("\n").split("\t")
            tmp = np.array([0, 0])
            for quiverId in collapse[line[1]]:
                tmp += ROIreads[quiverId]
            line = line+[str(tmp[0]), str(tmp[1])]
            out.append(line)
    with open(namespace.out, 'w') as File:
        for item in out:
            File.write("\t".join(item)+"\n")

'''
Descripttion: 
version: 
Author: zpliu
Date: 2021-02-15 15:40:24
LastEditors: zpliu
LastEditTime: 2021-02-16 08:52:21
@param: 
'''
import sys


def getConstitutiveIntronCoordinate(location):
    try:
        start, end = location.split(":")[2].split("-")
    except IndexError:
        start, end = location.split(":")[-1].split("-")
    strand = location.split(":")[-1]
    return int(start), int(end), strand


if __name__ == "__main__":
    ASkmerFile = sys.argv[1]
    conservedLocation = sys.argv[2]
    kmer = {}
    out = []
    with open(ASkmerFile, 'r') as File:
        for line in File:
            line = line.strip("\n").split("\t")
            if line[2] == '0':
                pass
            else:
                start1 = str(int(line[1].split(":")[-1].split("-")[0])+1)
                end1 = str(int(line[1].split(":")[-1].split("-")[1])-1)
                kmer[line[0]] = line[1].split(":")[0]+":"+start1+"-"+end1
    with open(conservedLocation, 'r') as File:
        for line in File:
            line = line.strip("\n").split("\t")
            start1, end1, strand1 = getConstitutiveIntronCoordinate(line[1])
            # no sense with starnd2
            start2, end2, strand2 = getConstitutiveIntronCoordinate(
                kmer[line[0]])
            tmp = sorted([start1, start2, end1, end2])
            if (tmp[2]-tmp[1])/(tmp[3]-tmp[0]) > 0.95:
                out.append("\t".join(line)+"\n")
            else:
                out.append(
                    line[0]+"\t"+line[1].split(":")[0]+":"+line[1].split(
                        ":")[1]+":"+str(start2)+"-"+str(end2)+":"+strand1+"\n"
                )
    with open(sys.argv[3], 'w') as File:
        for line in out:
            File.write(line)

'''
Descripttion: 
version: 
Author: zpliu
Date: 2021-01-04 10:28:37
LastEditors: zpliu
LastEditTime: 2021-01-04 16:20:48
@param: 
'''
import sys
from extractAScoordinate import getEventLength
import re


def getASeventBykmer(ASevent1, kmer, kmerLength, ASeventDict):
    '''
    @Descripttion: get the most identity AS event
    @param: genome1 AS event @str
    @param: genome2 kmer @str
    @param: kemer length @int
    @param: genome AS event @dict
    @return: 
    '''
    kmer_gene = kmer.split(";")[0].strip(">")
    kmerStart = int(kmer.split(":")[1].split("-")[0])
    kmerEnd = int(kmer.split(":")[1].split("-")[1])
    eventType = ASevent1.split(":")[0].split(";")[1]
    identity = 0
    conservedAS = ''
    try:
        ASevents = ASeventDict[kmer_gene]
    except KeyError:
      # homolog gene not contain AS
        return None, None
    for eventDict in ASevents:
        if re.search(eventType, eventDict['event']):
            if kmerStart > eventDict['coordinate'][1] or kmerEnd < eventDict['coordinate'][0]:
                # no intersect
                pass
            else:
                tmpcoordinate = sorted(
                    [kmerEnd, kmerStart, eventDict['coordinate'][1], eventDict['coordinate'][0]])

                tmp = tmpcoordinate[2]-tmpcoordinate[1]
                # the overlap lenth with at least 90% k-mer length
                if tmp/kmerLength > identity:
                    identity = tmp/kmerLength
                    conservedAS = eventDict['event']
        else:
            pass
    if identity < 0.9:
      # k-mer squence overlap with AS at leat 90%
        return None, None
    else:
        ASlength = getEventLength(conservedAS)['length']
        return conservedAS, ASlength


def getgeneEvent(ASFile):
    out = {}
    with open(ASFile, 'r') as File:
        for line in File:
            line = line.split("\t")
            out[line[1]] = out.get(line[1], [])
            out[line[1]].append(getEventLength(line[2]))
    return out


def getConstitutiveFeature(ConstitutiveFile):
    out = {}
    with open(ConstitutiveFile, 'r') as File:
        for line in File:
            line = line.strip("\n").split("\t")
            feature = {
                'length': int(line[2])-int(line[1])+1,
                'coordinate': [int(line[2]), int(line[1])],
                'event': line[3]+";"+line[-1]+":"+line[0]+":"+line[1]+"-"+line[2]+":"+line[4]
            }
            out[line[3]] = out.get(line[3], [])
            out[line[3]].append(getEventLength(line[2]))
    return out


def conservedAS(ASKmerFile, ASFile):
    out = []
    # all gene's AS event
    geneASDict = getgeneEvent(ASFile)
    with open(ASKmerFile, 'r') as File:
        for line in File:
            line = line.strip("\n").split("\t")
            if line[1]:
                ASevent, ASlength = getASeventBykmer(
                    line[0], line[1], int(line[3]), geneASDict)
                if ASevent:
                    # contain conserved AS event
                    out.append(
                        "\t".join([line[0], ASevent, str(ASlength), line[3]])+"\n")
                else:
                  # no conserved AS so used k-mer
                    out.append(
                        "\t".join([line[0], line[1], line[2], line[3]])+"\n")
            else:
                # no conserved kmer
                pass
    return out


if __name__ == "__main__":
    with open(sys.argv[3], 'w') as File:
        for item in conservedAS(sys.argv[1], sys.argv[2]):
            File.write(item)

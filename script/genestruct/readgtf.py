import numpy as np
import sys
import re
import operator


class gene(object):
    def __init__(self, geneName, stand, chromosome):
        self.geneName = geneName
        self.stand = stand
        self.transcripts = []
        self.chromosome = chromosome
        self.transcriptsObject = []

    def setTranscipts(self, transcriptId):
        if transcriptId not in self.transcripts:
            self.transcripts.append(transcriptId)
        else:
            return

    def setTranscriptsObject(self, transcriptDict):
        for item in self.transcripts:
            self.transcriptsObject.append(transcriptDict[item])

    def getMaxLengthTranscript(self):
        tmp = [i.getLength() for i in self.transcriptsObject]
        tmp2 = [i[1]-i[0]+1 for i in tmp]
        max_index, max_length = max(
            enumerate(tmp2), key=operator.itemgetter(1))  # 获取最大值及索引
        # [起始位置，终止位置]
        return tmp[max_index]

    def getTranscriptNum(self):
        # print(self.geneName)
        return len(self.transcriptsObject)

    def getTranscriptCoordinate(self, id):
        return self.transcriptsObject[id].getExonCoordinate()

    def getTranscriptName(self, id):
        return self.transcripts[id]

    def getChain(self):
        return self.stand

    def getAllspliceSite(self):
        tmp = []
        for item in self.transcriptsObject:
            tmp += item.getIntronCoordinate()  # 获取gene的所有intron坐标
        tmp2 = []
        [tmp2.append(item[0]) for item in tmp]
        [tmp2.append(item[1]) for item in tmp]
        return list(set(tmp2))

    def getAlltranscriptsObject(self):
        return self.transcriptsObject


class transcript:
    def __init__(self, transcriptName, stand):
        self.transcriptName = transcriptName
        self.stand = stand
        self.exons = []
        self.CDSs = []

    def setExcordinate(self, array, AnnotionType):
        if AnnotionType == 'exon':
            self.exons += array
        if AnnotionType == 'CDS':
            self.CDSs += array

    def getLength(self):
        if(len(self.exons) != 0):
            self.exons.sort()
            return [self.exons[0], self.exons[-1]]
        else:
            self.CDSs.sort()
            return [self.CDSs[0], self.CDSs[-1]]

    def getExonCoordinate(self):
        if(len(self.exons) != 0):
            return self.exons
        else:
            return self.CDSs

    def getIntronCoordinate(self):
        if(len(self.exons) == 2):
            return [[0, 0]]
        elif len(self.exons) > 2:
            self.exons.sort()
            tmp = []
            for i in range(1, len(self.exons)-1, 2):
                tmp.append([self.exons[i]+1, self.exons[i+1]-1])
            return tmp
        elif len(self.CDSs) == 2:
            return [[0, 0]]
        else:
            self.CDSs.sort()
            tmp = []
            for i in range(1, len(self.CDSs)-1, 2):
                tmp.append([self.CDSs[i]+1, self.CDSs[i+1]-1])
            return tmp

    def getTranscriptName(self):
        return self.transcriptName


def getGeneInfo(gtfFile):
    geneDict = {}
    transcriptDict = {}
    genePattern = r'gene_id "([^"]*)"'
    transcriptPattern = r'transcript_id "([^"]*)"'
    with open(gtfFile, 'r') as File:
        for line in File.readlines():
            line = line.strip("\n").split("\t")
            geneId = re.search(genePattern, line[8]).group(1)
            transcriptId = re.search(transcriptPattern, line[8]).group(1)
            try:
                transcriptDict[transcriptId].setExcordinate(
                    [int(line[3]), int(line[4])], line[2])
            except:
                transcriptDict[transcriptId] = transcript(
                    transcriptId, line[6])
                transcriptDict[transcriptId].setExcordinate(
                    [int(line[3]), int(line[4])], line[2])  # 存坐标轴信息
            try:
                geneDict[geneId].setTranscipts(transcriptId)
            except KeyError:
                geneDict[geneId] = gene(geneId, line[6], line[0])
                geneDict[geneId].setTranscipts(transcriptId)  # 存转录本信息
    # 将转录本对象存进gene对象中
    for item in geneDict.keys():
        geneDict[item].setTranscriptsObject(transcriptDict)
    return geneDict


if __name__ == "__main__":
    geneDict = getGeneInfo(sys.argv[1])
    for key in geneDict.keys():
        print(geneDict[key].getTranscriptCoordinate(0))
        break

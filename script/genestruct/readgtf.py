# coding:utf-8
import numpy as np
import sys
import re
import copy
from readFastaFromBed import reverseSequence


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
    def __init__(self, transcriptName, stand, chromosome, geneId):
        self.transcriptName = transcriptName
        self.stand = stand
        self.chromosome = chromosome
        self.geneId = geneId
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

    def getgeneName(self):
        return self.geneId

    def getTranscriptName(self):
        return self.transcriptName

    def getcDNAlength(self):
        if(len(self.exons) != 0):
            self.exons.sort()
            tmp = 0
            for i in range(0, len(self.exons), 2):
                tmp += self.exons[i+1]-self.exons[i]+1
            return tmp
        else:
            self.CDSs.sort()
            tmp = 0
            for i in range(0, len(self.CDSs), 2):
                tmp += self.CDSs[i+1]-self.CDSs[i]+1
            return tmp

    def getcDNAlengthBycoordinate(self, start, end, refStart, refEnd, PacBioexon):
        # 根据起始密码子和终止密码子的位置提取参考基因组对应转录本的cDNA序列长度
        # 根据终止密码子和起始密码子可以获得对应参考基因的cDNA序列，看是否是3的倍数
        # 这里需要使用到深拷贝 ,直接赋值会改变原有的对象
        def caculateLength(start, end, exons):  # 包括end这个坐标
            length = 0
            exons.sort()
            for i in range(0, len(exons), 2):  # 计算一定范围的cDNA碱基数
                if exons[i+1] < start or exons[i] > end:
                    continue
                elif exons[i] < start and exons[i+1] <= end:
                    length += exons[i+1]-start+1
                elif exons[i+1] > end and exons[i] >= start:
                    length += end-exons[i]+1
                elif exons[i] >= start and exons[i+1] <= end:
                    length += exons[i+1]-exons[i]+1
                else:
                    length += end-start+1
            return length
        if(len(self.exons) != 0):
            if start == refStart and end > refEnd and self.stand == "+":
                return caculateLength(start, refEnd, PacBioexon)
            if start == refStart and end > refEnd and self.stand == "-":
                return caculateLength(end, start, self.exons)
            if start == refStart and end < refEnd and self.stand == "+":
                return caculateLength(start, end, self.exons)
            if start == refStart and end < refEnd and self.stand == "-":
                return caculateLength(refEnd, start, PacBioexon)
            if end == refEnd and start > refStart and self.stand == "+":
                return caculateLength(refStart, start-1, self.exons)
            if end == refEnd and start > refStart and self.stand == "-":
                return caculateLength(refStart+1, start, PacBioexon)
            if end == refEnd and start < refStart and self.stand == "+":
                return caculateLength(start, refStart-1, PacBioexon)
            if end == refEnd and start < refStart and self.stand == "-":
                return caculateLength(start+1, refStart, self.exons)
        else:
            if start == refStart and end > refEnd and self.stand == "+":
                return caculateLength(start, refEnd, PacBioexon)
            if start == refStart and end > refEnd and self.stand == "-":
                return caculateLength(end, start, self.CDSs)
            if start == refStart and end < refEnd and self.stand == "+":
                return caculateLength(start, end, self.CDSs)
            if start == refStart and end < refEnd and self.stand == "-":
                return caculateLength(refEnd, start, PacBioexon)
            if end == refEnd and start > refStart and self.stand == "+":
                return caculateLength(refStart, start-1, self.CDSs)
            if end == refEnd and start > refStart and self.stand == "-":
                return caculateLength(refStart+1, start, PacBioexon)
            if end == refEnd and start < refStart and self.stand == "+":
                return caculateLength(start, refStart-1, PacBioexon)
            if end == refEnd and start < refStart and self.stand == "-":
                return caculateLength(start+1, refStart, self.CDSs)

    def getExonCoordinate(self):
        if(len(self.exons) != 0):
            return self.exons
        else:
            return self.CDSs

    def getIntronCoordinate(self):
        '''
        根据外显子的注释信息，获取对应的内含子坐标
        '''
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

    def getAbsoluteCoordinate(self, start, end):
        '''
        根据转录本ORF序列的相对位置，定位到基因组上的绝对坐标
        + 如果gene位于正链

        + 如果gene位于负链
        '''
        if self.stand == "+":
            if len(self.exons) != 0:
                self.exons.sort()
                for i in range(0, len(self.exons), 2):
                    if start <= self.exons[i+1]-self.exons[i]+1:
                        AbsoulateStart = self.exons[i]+start-1
                        break
                    else:
                        start = start - \
                            (self.exons[i+1]-self.exons[i]+1)  # 获取绝对起始位置
                isofomrlength = 0
                for i in range(0, len(self.exons), 2):
                    isofomrlength += self.exons[i+1]-self.exons[i]+1
                    if isofomrlength < end:
                        continue
                    else:
                        AbsoulateEnd = self.exons[i+1]-(isofomrlength-end)
                        break
            else:
                self.CDSs.sort()
                for i in range(0, len(self.CDSs), 2):
                    if start <= self.CDSs[i+1]-self.CDSs[i]+1:
                        AbsoulateStart = self.CDSs[i]+start-1
                        break
                    else:
                        start = start - \
                            (self.CDSs[i+1]-self.CDSs[i]+1)  # 获取绝对起始位置
                isofomrlength = 0
                for i in range(0, len(self.CDSs), 2):
                    isofomrlength += self.CDSs[i+1]-self.CDSs[i]+1
                    if isofomrlength < end:
                        continue
                    else:
                        AbsoulateEnd = self.CDSs[i+1]-(isofomrlength-end)
                        break

        else:
            if len(self.exons) != 0:
                self.exons.sort()
                for i in range(len(self.exons)-1, 0, -2):
                    if self.exons[i]-self.exons[i-1]+1 > start:
                        AbsoulateStart = self.exons[i]-start+1
                        break  # 找到了就得停下来，不然会一直找下去
                    else:
                        start = start-(self.exons[i]-self.exons[i-1]+1)
                isofomrlength = 0
                for i in range(len(self.exons)-1, 0, -2):
                    isofomrlength += self.exons[i]-self.exons[i-1]+1
                    if isofomrlength < end:
                        continue
                    else:
                        AbsoulateEnd = self.exons[i-1]+(isofomrlength-end)
                        break
            else:
                self.CDSs.sort()
                for i in range(len(self.CDSs)-1, 0, -2):
                    if self.CDSs[i]-self.CDSs[i-1]+1 > start:
                        AbsoulateStart = self.CDSs[i]-start+1
                        break
                    else:
                        start = start-(self.CDSs[i]-self.CDSs[i-1]+1)
                isofomrlength = 0
                for i in range(len(self.CDSs)-1, 0, -2):
                    isofomrlength += self.CDSs[i]-self.CDSs[i-1]+1
                    if isofomrlength < end:
                        continue
                    else:
                        AbsoulateEnd = self.CDSs[i-1]+(isofomrlength-end)
                        break
        return [AbsoulateStart, AbsoulateEnd]

    def getcDNAsequence(self, genomeSequence):
        '''
        根据exon的注释信息获取对应的cDNA序列
        '''
        sequence = ''
        if(len(self.exons) != 0):
            self.exons.sort()
            for i in range(0, len(self.exons), 2):
                sequence += genomeSequence[self.chromosome][self.exons[i] -
                                                            1:self.exons[i+1]]
        else:
            self.CDSs.sort()
            for i in range(1, len(self.CDSs), 2):
                sequence += genomeSequence[self.chromosome][self.CDSs[i] -
                                                            1:self.CDSs[i+1]]
        if self.stand == "+":
            return ">"+self.transcriptName+"\n"+sequence
        else:
            return ">"+self.transcriptName+"\n"+reverseSequence(sequence)


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
                    transcriptId, line[6], line[0], geneId)
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


def getTranscriptInfo(gtfFile):
    '''
    从gtf文件中获取每条转录本的注释信息
    '''
    transcriptDict = {}
    genePattern = r'gene_id "([^"]*)"'
    transcriptPattern = r'transcript_id "([^"]*)"'
    with open(gtfFile, 'r') as File:
        for line in File.readlines():
            line = line.strip("\n").split("\t")
            geneId = re.search(genePattern, line[8]).group(1)
            transcriptId = re.search(transcriptPattern, line[8]).group(1)
            # print(transcriptId, end="\t")
            # print(geneId, end="\t")
            try:
                transcriptDict[transcriptId].setExcordinate(
                    [int(line[3]), int(line[4])], line[2])
            except:
                transcriptDict[transcriptId] = transcript(
                    transcriptId, line[6], line[0], geneId)
                transcriptDict[transcriptId].setExcordinate(
                    [int(line[3]), int(line[4])], line[2])  # 存坐标轴信息
    return transcriptDict


if __name__ == "__main__":
    geneDict = getGeneInfo(sys.argv[1])
    for key in geneDict.keys():
        print(geneDict[key].getTranscriptCoordinate(0))
        break

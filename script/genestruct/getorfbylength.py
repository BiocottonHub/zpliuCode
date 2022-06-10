'''
Descripttion: 
version: 
Author: zpliu
Date: 2020-12-28 19:25:14
LastEditors: zpliu
LastEditTime: 2020-12-28 19:28:36
@param: isoform Annotion gtf file @str
@param: EMBOSS all predict CDS sequence file @str
@paran: out file of CDS sequence @str
'''
# coding:utf-8


import sys
import re
from readgtf import getTranscriptInfo


class isoformORF:
    def __init__(self, isoformName, geneName):
        self.isoformName = isoformName
        self.geneName = geneName
        self.ORFCoordinate = []
        self.ORFAbsoluteCoordinate = []
        self.allORF = []

    def getgeneName(self):
        return self.geneName

    def setCoordinate(self, start, end):
        self.ORFCoordinate.append([start, end])

    def setORF(self, sequence):
        self.allORF.append(sequence)

    def getORF(self):
        # 获取最长的那个ORF作为转录本的ORF信息
        tmp = [i[1]-i[0] for i in self.ORFCoordinate]
        index = tmp.index(max(tmp))
        # print(len(self.allORF), end="\t")
        # print(len(self.ORFCoordinate), end="\t")
        # print(self.isoformName)
        sequence = self.allORF[index]
        return {'name': self.isoformName, 'Relativestart': self.ORFCoordinate[index][0], 'Relativeend': self.ORFCoordinate[index][1], 'sequence': sequence}

    def getAbsoulteORF(self):
        return self.ORFAbsoluteCoordinate

    def setAbsoulteORF(self, start, end):
        self.ORFAbsoluteCoordinate.append(start)
        self.ORFAbsoluteCoordinate.append(end)


def getisoformORF(gtfFile, ORFFile):
    '''
    根据转录本的注释信息和ORF的预测信息
    获取起始密码子和终止密码的绝对坐标和相对坐标
    '''
    isoformDict = getTranscriptInfo(gtfFile)  # 获取转录本的坐标信息
    isoformORFDict = {}  # 存储转录本的ORF注释信息
    with open(ORFFile, 'r') as File:
        sequence = ''
        for line in File.readlines():
            if re.match('^>', line):
                try:
                    isoformORFDict[isoformName].setORF(sequence)  #
                    sequence = ''
                except NameError:  # 第一行用于跳过
                    pass  # 啥也不用干
                isoformName = re.split(r'_[0-9]+', line)[0].strip(">")
                start = int(
                    re.search('\[([^\]]*)\]', line).group(1).split("-")[0].strip(" "))
                end = int(
                    re.search('\[([^\]]*)\]', line).group(1).split("-")[1].strip(" "))
                try:
                    isoformORFDict[isoformName].setCoordinate(start, end)
                except KeyError:
                    # 获取转录本对应的gene编号
                    geneName = isoformDict[isoformName].getgeneName()
                    isoformORFDict[isoformName] = isoformORF(
                        isoformName, geneName)
                    isoformORFDict[isoformName].setCoordinate(start, end)
            else:
                sequence += line.strip("\n")
        isoformORFDict[isoformName].setORF(sequence)  # 最后一行的序列没有加进去
    '''
    根据isoform的注释信息；更新ORF的绝对坐标信息
    '''
    for key in isoformORFDict.keys():
        relativeStart = isoformORFDict[key].getORF()['Relativestart']
        relativeEnd = isoformORFDict[key].getORF()['Relativeend']
        # print(key, end="\t")
        # print(relativeStart, end="\t")
        # print(relativeEnd)
        [AbsoluteStart, AbsoluteEnd] = isoformDict[key].getAbsoluteCoordinate(
            relativeStart, relativeEnd)
        isoformORFDict[key].setAbsoulteORF(
            AbsoluteStart, AbsoluteEnd)  # 更新了ORF的绝对坐标信息
    return isoformORFDict


if __name__ == "__main__":
    isoformDict = getisoformORF(sys.argv[1], sys.argv[2])
    with open(sys.argv[3], 'w') as File:
        for key in isoformDict.keys():
            # print(isoformDict[key].getORF())
            # print(isoformDict[key].getAbsoulteORF())
            out = isoformDict[key].getORF()
            File.write(">"+out['name']+"\n"+out['sequence']+"\n")

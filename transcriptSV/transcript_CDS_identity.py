'''
Descripttion: 
version: 
Author: zpliu
Date: 2020-12-23 22:28:57
LastEditors: zpliu
LastEditTime: 2020-12-24 09:48:25
@param: 
'''
import argparse
from muscle import proteinMuscle
from getProductiveTranscript import getProdectiveTranscript
from fun.readFasta import readFastaFile
from itertools import combinations


class IsoformMessage(object):
    ''' 
    @Descripttion: the PacBio Isoform Message
    @param: isoform name @str
    @param: CDS sequence @str
    '''
    #__slots__ = ('isoform', 'sequence', 'sequenceLen')
    ##############################
    #isoform @str
    #############################
    @property
    def isoform(self):
        return self.__isoform

    @isoform.setter
    def isoform(self,  isoformName: str):
        self.__isoform = isoformName
    ##############################
    #sequence @str
    #sequenceLen @int
    #############################
    @property
    def sequence(self):
        return self.__sequence

    @sequence.setter
    def sequence(self, fastaSequence):
        self.__sequence = ">{name}\n{fasta}".format(
            name=self.__isoform, fasta=fastaSequence)
        self.__sequenceLen = len(fastaSequence)

    @property
    def sequenceLen(self):
        return self.__sequenceLen


def getGeneMessage(geneIsoformFile, IsoformSequenceFile, prefix):
    '''
    @Descripttion: According the gene and Isoform message to designer the dict
    @param: gene and isoform @str
    @param: isoform sequence @str
    @param: isoform's prefix @str
    @return: gene and isoform Dict @dict{@dict}
    '''
    IsoformSequenceDict = readFastaFile(IsoformSequenceFile)
    out = {}
    with open(geneIsoformFile, 'r') as File:
        for line in File:
            line = line.strip("\n").split("\t")
            isoformName = "{prefix}^{name}".format(
                prefix=prefix, name=line[1])
            out[line[0]] = out.get(line[0], {})
            out[line[0]][isoformName] = IsoformMessage()
            out[line[0]][isoformName].isoform = isoformName
            out[line[0]][isoformName].sequence = IsoformSequenceDict[line[1]]
    return out


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-homolog', help='homolog gene File')
    parser.add_argument('-fasta1', help='A transcripts cDNA')
    parser.add_argument('-fasta2', help='B transcripts cDNA')
    parser.add_argument('-RNAseq1', help='A transcripts RNA-seq')
    parser.add_argument('-RNAseq2', help='B transcripts RNA-seq')
    parser.add_argument('-prex1', help='A transcripts prefix')
    parser.add_argument('-prex2', help='B transcripts prefix')
    parser.add_argument('-out', help='out file')
    nameSpace = parser.parse_args()
    AgeneIsoform = getGeneMessage(
        nameSpace.RNAseq1, nameSpace.fasta1, nameSpace.prex1)
    BgeneIsoform = getGeneMessage(
        nameSpace.RNAseq2, nameSpace.fasta2, nameSpace.prex2)
    out = []
    with open(nameSpace.homolog, 'r') as File:
        for line in File:
            line = line.strip("\n").split("\t")
            AgeneId = line[2]  # A gene column
            BgeneId = line[3]  # B gene column
            if AgeneId in AgeneIsoform and BgeneId in BgeneIsoform:
                AisoformClassDict = AgeneIsoform[AgeneId]
                BisoformClassDict = BgeneIsoform[BgeneId]
                AllIsoform = []
                [AllIsoform.append(i) for i in AisoformClassDict.values()]
                [AllIsoform.append(i) for i in BisoformClassDict.values()]
                for item in combinations(AllIsoform, 2):
                  # to get full combnation of isoforms
                  # filter with sequenc is same
                    if item[0].sequenceLen == item[1].sequenceLen:
                        idnetity = str(proteinMuscle(
                            item[0].sequence, item[1].sequence)/item[0].sequenceLen)
                        out.append(
                            "{}\t{}\t{}\t{}\t{}\t{}\n".format(
                                AgeneId, BgeneId, item[0].isoform, item[1].isoform, item[0].sequenceLen, idnetity
                            )
                        )
                    else:
                        pass
    with open(nameSpace.out, 'w') as File:
        for item in out:
            File.write(item)

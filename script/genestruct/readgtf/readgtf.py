from readFasta.readFastaFromBed import reverseSequence
import re

pattern = 'transcript_id \"([^\"]*)\"'  # Match isform id from gtf file


class isformClass:
    def __init__(self, isformId, chain, chrosome):
        self.isformId = isformId,
        self.chromosome = chrosome
        self.exonCoordinate = []
        self.chain = chain
        self.CDS = []

    def setExon(self, coordinate1, coordinate2):
        self.exonCoordinate.append(coordinate1)
        self.exonCoordinate.append(coordinate2)

    def setCDS(self, coordinate1, coordinate2):
        self.CDS.append(coordinate1)
        self.CDS.append(coordinate2)

    def getCDNAsequence(self, genomeSequence):

        sequence = ''
        if(len(self.exonCoordinate) != 0):
            self.exonCoordinate.sort()
            # print(self.exonCoordinate, end="\t")
            # print(self.isformId)
            for i in range(1, len(self.exonCoordinate), 2):
                sequence += genomeSequence[self.chromosome][self.exonCoordinate[i-1] -
                                                            1:self.exonCoordinate[i]]
        else:
            self.CDS.sort()
            for i in range(1, len(self.CDS), 2):
                sequence += genomeSequence[self.chromosome][self.CDS[i-1] -
                                                            1:self.CDS[i]]
        if self.chain == "+":
            return sequence
        else:
            return reverseSequence(sequence)

    def getExoncoordinate(self):

        if(len(self.exonCoordinate) != 0):
            self.exonCoordinate.sort()
            return [self.exonCoordinate, self.chain, self.chromosome]
        else:
            self.CDS.sort()
            return [self.CDS, self.chain, self.chromosome]

    def finderFlankExonSequencTag(self, ASEventId, EventStart, EventEnd, genomeSequence):
        '''
        According AS coordinate find flanking exon sequence Tag
        '''
        # self.printExonCoordinate()
        if len(self.exonCoordinate) != 0:
            junctionIndex1 = 0
            junctionIndex2 = len(self.exonCoordinate)-1
            while(self.exonCoordinate[junctionIndex1] < EventStart and junctionIndex1 < junctionIndex2):
                junctionIndex1 = junctionIndex1+1
            while(self.exonCoordinate[junctionIndex2] > EventEnd and junctionIndex2 > 0):
                junctionIndex2 = junctionIndex2-1
            if(junctionIndex1 == len(self.exonCoordinate)-1 or junctionIndex2 == 0):
                print("超出isfrom范围事件\t"+ASEventId)
                return ""
            FEST5Sequence = ''
            FEST3Sequence = ''
            try:
                if self.chain == "+" and self.exonCoordinate[junctionIndex2] - self.exonCoordinate[junctionIndex2-1] >= 299:
                    FEST5Sequence = genomeSequence[self.chromosome][(
                        self.exonCoordinate[junctionIndex2] - 300):self.exonCoordinate[junctionIndex2]]
                if self.chain == "+" and self.exonCoordinate[junctionIndex2] - self.exonCoordinate[junctionIndex2-1] < 299:
                    FEST5Sequence = genomeSequence[self.chromosome][self.exonCoordinate[junctionIndex2-1] -
                                                                    1:self.exonCoordinate[junctionIndex2]]
                if self.chain == "+" and self.exonCoordinate[junctionIndex1+1] - self.exonCoordinate[junctionIndex1] >= 299:
                    FEST3Sequence = genomeSequence[self.chromosome][(
                        self.exonCoordinate[junctionIndex1] - 1):(self.exonCoordinate[junctionIndex1]+299)]
                if self.chain == "+" and self.exonCoordinate[junctionIndex1+1] - self.exonCoordinate[junctionIndex1] < 299:
                    FEST3Sequence = genomeSequence[self.chromosome][(
                        self.exonCoordinate[junctionIndex1] - 1):(self.exonCoordinate[junctionIndex1+1])]
                if self.chain == "-" and self.exonCoordinate[junctionIndex1+1] - self.exonCoordinate[junctionIndex1] >= 299:
                    FEST5Sequence = reverseSequence(
                        genomeSequence[self.chromosome][(self.exonCoordinate[junctionIndex1]-1):(self.exonCoordinate[junctionIndex1]+299)])
                if self.chain == "-" and self.exonCoordinate[junctionIndex1+1] - self.exonCoordinate[junctionIndex1] < 299:
                    FEST5Sequence = reverseSequence(
                        genomeSequence[self.chromosome][(self.exonCoordinate[junctionIndex1]-1):self.exonCoordinate[junctionIndex1+1]])
                if self.chain == "-" and self.exonCoordinate[junctionIndex2] - self.exonCoordinate[junctionIndex2-1] >= 299:
                    FEST3Sequence = reverseSequence(
                        genomeSequence[self.chromosome][(self.exonCoordinate[junctionIndex2]-300):self.exonCoordinate[junctionIndex2]])
                if self.chain == "-" and self.exonCoordinate[junctionIndex2] - self.exonCoordinate[junctionIndex2-1] < 299:
                    FEST3Sequence = reverseSequence(genomeSequence[self.chromosome][(
                        self.exonCoordinate[junctionIndex2-1]-1):self.exonCoordinate[junctionIndex2]])
            except IndexError as reason:
                print(str(self.isformId)+"AS坐标有误"+str(reason))
            except KeyError as reason:
                print(str(self.chromosome)+"对应染色体在基因组文件中不存在"+str(reason))
            ASEventId = ASEventId+"-"+self.chain+"-"
            if(len(FEST5Sequence) >= 20 and len(FEST3Sequence) >= 20):
                return(">"+ASEventId+"UpFEST\n"+FEST5Sequence+"\n"+">"+ASEventId+"DownFEST\n"+FEST3Sequence+"\n")
            else:
                return ''
        else:
            junctionIndex1 = 0
            junctionIndex2 = len(self.CDS)-1
            while(self.CDS[junctionIndex1] < EventStart and junctionIndex1 < junctionIndex2):
                junctionIndex1 = junctionIndex1+1
            while(self.CDS[junctionIndex2] > EventEnd and junctionIndex2 > 0):
                junctionIndex2 = junctionIndex2-1
            if(junctionIndex1 == len(self.CDS)-1 or junctionIndex2 == 0):
                print("超出isfrom范围事件\t"+ASEventId)
                return ""
            FEST5Sequence = ''
            FEST3Sequence = ''
            try:
                if self.chain == "+" and self.CDS[junctionIndex2] - self.CDS[junctionIndex2-1] >= 299:
                    FEST5Sequence = genomeSequence[self.chromosome][(
                        self.CDS[junctionIndex2] - 300):self.CDS[junctionIndex2]]
                if self.chain == "+" and self.CDS[junctionIndex2] - self.CDS[junctionIndex2-1] < 299:
                    FEST5Sequence = genomeSequence[self.chromosome][self.CDS[junctionIndex2-1] -
                                                                    1:self.CDS[junctionIndex2]]
                if self.chain == "+" and self.CDS[junctionIndex1+1] - self.CDS[junctionIndex1] >= 299:
                    FEST3Sequence = genomeSequence[self.chromosome][(
                        self.CDS[junctionIndex1] - 1):(self.CDS[junctionIndex1]+299)]
                if self.chain == "+" and self.CDS[junctionIndex1+1] - self.CDS[junctionIndex1] < 299:
                    FEST3Sequence = genomeSequence[self.chromosome][(
                        self.CDS[junctionIndex1] - 1):(self.CDS[junctionIndex1+1])]
                if self.chain == "-" and self.CDS[junctionIndex1+1] - self.CDS[junctionIndex1] >= 299:
                    FEST5Sequence = reverseSequence(
                        genomeSequence[self.chromosome][(self.CDS[junctionIndex1]-1):(self.CDS[junctionIndex1]+299)])
                if self.chain == "-" and self.CDS[junctionIndex1+1] - self.CDS[junctionIndex1] < 299:
                    FEST5Sequence = reverseSequence(
                        genomeSequence[self.chromosome][(self.CDS[junctionIndex1]-1):self.CDS[junctionIndex1+1]])
                if self.chain == "-" and self.CDS[junctionIndex2] - self.CDS[junctionIndex2-1] >= 299:
                    FEST3Sequence = reverseSequence(
                        genomeSequence[self.chromosome][(self.CDS[junctionIndex2]-300):self.CDS[junctionIndex2]])
                if self.chain == "-" and self.CDS[junctionIndex2] - self.CDS[junctionIndex2-1] < 299:
                    FEST3Sequence = reverseSequence(genomeSequence[self.chromosome][(
                        self.CDS[junctionIndex2-1]-1):self.CDS[junctionIndex2]])
            except IndexError as reason:
                print(str(self.isformId)+"AS坐标有误"+str(reason))
            except KeyError as reason:
                print(str(self.chromosome)+"对应染色体在基因组文件中不存在"+str(reason))
            ASEventId = ASEventId+"-"+self.chain+"-"
            if(len(FEST5Sequence) >= 20 and len(FEST3Sequence) >= 20):
                return(">"+ASEventId+"UpFEST\n"+FEST5Sequence+"\n"+">"+ASEventId+"DownFEST\n"+FEST3Sequence+"\n")
            else:
                return ''

    def isoformAS(self, start, end):
        '''
        判断哪个一个转录本发生AS,外显子坐标与事件坐标没有交集，则没有发生AS
        True: 该转录本发生AS
        False:该转录本没有发生AS
        '''
        if self.exonCoordinate == []:
            for index in range(0, len(self.CDS), 2):
                if(self.CDS[index] > end or self.CDS[index+1] < start):
                    continue
                else:
                    return True
        else:
            for index in range(0, len(self.exonCoordinate), 2):
                if(self.exonCoordinate[index] > end or self.exonCoordinate[index+1] < start):
                    continue
                else:
                    return True
        return False

    def getASPosition(self, start, end):
        tmp = []
        if self.exonCoordinate == []:
            for index in range(0, len(self.CDS), 2):
                if(self.CDS[index] > end or self.CDS[index+1] < start):
                    continue
                else:
                    tmp.append([self.CDS[index], self.CDS[index+1]])

        else:
            for index in range(0, len(self.exonCoordinate), 2):
                if(self.exonCoordinate[index] > end or self.exonCoordinate[index+1] < start):
                    continue
                else:
                    tmp.append([self.exonCoordinate[index],
                                self.exonCoordinate[index+1]])
        return tmp

    def getIntronPosition(self):
        self.exonCoordinate.sort()
        self.CDS.sort()
        tmp = []
        if self.exonCoordinate == []:
            for index in range(1, len(self.CDS)-1, 2):

                tmp.append(
                    [self.chromosome, self.CDS[index]+1, self.CDS[index+1]-1])

        else:
            for index in range(1, len(self.exonCoordinate)-1, 2):
                tmp.append([self.chromosome, self.exonCoordinate[index]+1,
                            self.exonCoordinate[index+1]-1])
        return tmp


def readIsformFromgtf(fileName):
    isoforms = {}
    with open(fileName, 'r') as gtfFile:
        lines = gtfFile.readlines()
    for item in lines:
        item = item.split("\t")
        isformId = re.search(pattern, item[8]).group(1)
        if isformId not in isoforms and re.match('exon', item[2]):
            isoforms[isformId] = isformClass(isformId, item[6], item[0])
            isoforms[isformId].setExon(int(item[3]), int(item[4]))
        elif isformId not in isoforms and re.match('CDS', item[2]):
            isoforms[isformId] = isformClass(isformId, item[6], item[0])
            isoforms[isformId].setCDS(int(item[3]), int(item[4]))
        elif re.match('exon', item[2]):
            isoforms[isformId].setExon(int(item[3]), int(item[4]))
        elif re.match('CDS', item[2]):
            isoforms[isformId].setCDS(int(item[3]), int(item[4]))
    return isoforms

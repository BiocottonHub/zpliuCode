'''
Descripttion: 
version: 
Author: zpliu
Date: 2021-03-11 11:24:24
LastEditors: zpliu
LastEditTime: 2021-03-14 17:50:21
#! @param: read gtf file from gffcompare
'''
import sys
import pickle


class transcript(object):
    def __init__(self, AttributesDict: dict):
        super().__init__()
        # * initialization the attributions
        '''
        @param: exon@list
        @param: transcript_id@str
        @param: cmp_ref@str
        @param: ref_gene_id@str
        @param: cmp_ref_gene@str
        @param: stand@str
        @param: start@int
        @param: end@str
        @param: xloc@str
        @param: class_code@str
        @param: gene_id@str
        '''
        self.exon = AttributesDict['exon']
        self.transcript_id = AttributesDict['transcript_id']
        self.cmp_ref = AttributesDict['cmp_ref']
        self.ref_gene_id = AttributesDict['ref_gene_id']
        self.cmp_ref_gene = AttributesDict['cmp_ref_gene']
        self.stand = AttributesDict['stand']
        self.start = AttributesDict['start']
        self.end = AttributesDict['end']
        self.xloc = AttributesDict['xloc']
        self.class_code = AttributesDict['class_code']
        self.gene_id = AttributesDict['gene_id']
        self.chromsome = AttributesDict['chromsome']

    def set_exonCoordinate(self, coordinate: tuple):
        # * set exon coordinate
        self.exon.append(coordinate[0])
        self.exon.append(coordinate[1])

    def get_GeneId(self):
        return self.gene_id

    def get_intronCoordinate(self):
        sortedExon = sorted(self.exon)
        intronCoordinate = []
        for index in range(2, len(sortedExon)-1, 2):
          # * intron tubple with start and end
            intronCoordinate.append(
                (self.chromsome,
                 sortedExon[index-1]+1,
                 sortedExon[index]-1
                 )
            )
        return intronCoordinate

    def get_transcriptGtf(self):
        # * get the message of transcript line in gtf file
        line0 = self.chromsome
        line1 = 'StringTie'
        line2 = 'transcript'
        line3 = str(self.start)
        line4 = str(self.end)
        line5 = '.'
        line6 = self.stand
        line7 = '.'
        tmp1 = 'transcript_id \"{}\"; gene_id \"{}\"; xloc \"{}\"; class_code \"{}\";'.format(
            self.transcript_id,
            self.gene_id,
            self.xloc,
            self.class_code
        )
        if self.cmp_ref_gene:
            tmp2 = 'cmp_ref_gene \"{}\";'.format(self.cmp_ref_gene)
        else:
            tmp2 = 'cmp_ref_gene \"{}\";'.format(self.ref_gene_id)
        return "\t".join([line0, line1, line2, line3, line4, line5, line6, line7, tmp1+tmp2])+"\n"

    def get_exonGtf(self):
        out = ''
        line0 = self.chromsome
        line1 = 'StringTie'
        line2 = 'exon'
        line5 = '.'
        line6 = self.stand
        line7 = '.'
        for index in range(0, len(self.exon), 2):
            line3 = str(self.exon[index])
            line4 = str(self.exon[index+1])
            tmp = 'transcript_id \"{}\"; gene_id \"{}\"; exon_number \"{}\";'.format(
                self.transcript_id,
                self.gene_id,
                str(int((index+2)/2))
            )
            out += "\t".join([line0, line1, line2, line3,
                              line4, line5, line6, line7, tmp])+"\n"
        return out

    def __str__(self):
        ''' 
        # get description of class object
        '''
        attributeStr = ",".join("{}={}".format(k, getattr(self, k))
                                for k in self.__dict__.keys())
        return "[class name:{}\n{}]".format(self.__class__.__name__, attributeStr)


def getAttributesDict(attributesStr: str):
    # * get attributions of transcripts
    # todo
    #  gene_name "Ghir_A01G000020"; xloc "XLOC_000003"; ref_gene_id
    # "Ghir_A01G000020"; cmp_ref "Ghir_A01G000020.1";
    # class_code "o"; tss_id "TSS6"
    # todo
    return [i.strip(" ").split(" ") for i in
            attributesStr.replace('"', "").split(";")[:-1]]


def readgtfFile(gtfFile):
    transcriptList = {}
    with open(gtfFile, 'r') as File:
        for line in File:
            if line.startswith('#'):
                pass
            else:
                line = line.strip("\n").split("\t")
                if line[2] == "transcript":
                    # * set attributes of transcript
                    attributesDict = dict(getAttributesDict(line[8]))
                    attributesDict['exon'] = []
                    attributesDict['start'] = int(line[3])
                    attributesDict['end'] = int(line[4])
                    attributesDict['stand'] = line[6]
                    attributesDict['chromsome'] = line[0]
                    if "cmp_ref" not in attributesDict:
                      # ! transcript in intergenic
                        continue
                    # ? different compared to reference gene
                    attributesDict['ref_gene_id'] = attributesDict.get(
                        'ref_gene_id', None)
                    attributesDict['cmp_ref_gene'] = attributesDict.get(
                        'cmp_ref_gene', None)
                    transcriptList[attributesDict['transcript_id']] = transcriptList.get(
                        attributesDict['transcript_id'], transcript(attributesDict))
                else:
                    try:
                        attributesList = dict(getAttributesDict(line[8]))
                        transcriptList[attributesDict['transcript_id']].set_exonCoordinate(
                            (int(line[3]), int(line[4]))
                        )
                    except KeyError:
                      # ! transcript in intergenic
                        continue
    return transcriptList


if __name__ == "__main__":
    gffFile = sys.argv[1]
    # readgtfFile(gffFile)
    with open(sys.argv[2], 'wb') as File:
        # for key in readgtfFile(gffFile).keys():
        #     File.write(key+"\n")
        pickle.dump(readgtfFile(gffFile), File)

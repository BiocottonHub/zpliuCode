'''
Descripttion: 
version: 
Author: zpliu
Date: 2020-12-22 20:14:56
LastEditors: zpliu
LastEditTime: 2020-12-22 20:44:44
@param: 
'''
import argparse
from muscle import runMuscle
from getProductiveTranscript import getProdectiveTranscript
from fun.readFasta import readFastaFile


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-homolog', help='homolog gene File')
    parser.add_argument('-fasta1', help='A transcripts cDNA')
    parser.add_argument('-fasta2', help='B transcripts cDNA')
    parser.add_argument('-RNAseq1', help='A transcripts RNA-seq')
    parser.add_argument('-RNAseq2', help='B transcripts RNA-seq')
    parser.add_argument('-out', help='out file')
    nampeSpace = parser.parse_args()
    ProductiveIsoform1 = getProdectiveTranscript(nampeSpace.RNAseq1)
    ProductiveIsoform2 = getProdectiveTranscript(nampeSpace.RNAseq2)
    isoform1Fasta = readFastaFile(nampeSpace.fasta1)
    isoform2Fasta = readFastaFile(nampeSpace.fasta2)
    out = []
    with open(nampeSpace.homolog, 'r') as File:
        for line in File:
            line = line.strip("\n").split("\t")
            AgeneId = line[2]  # A gene column
            BgeneId = line[3]  # B gene column
            if AgeneId in ProductiveIsoform1 and BgeneId in ProductiveIsoform2:
                AtranscriptMessage = ProductiveIsoform1[AgeneId]
                BtranscriptMessage = ProductiveIsoform2[BgeneId]
                Afasta = ">"+AtranscriptMessage[0]+"\n" + \
                    isoform1Fasta[AtranscriptMessage[0]]+"\n"
                Bfasta = ">"+BtranscriptMessage[0]+"\n" + \
                    isoform2Fasta[BtranscriptMessage[0]]+"\n"
                try:
                    delteCode = runMuscle(Afasta, Bfasta)
                except:
                    print(line)
                tmp = "\t".join(
                    [AgeneId, BgeneId, AtranscriptMessage[0], BtranscriptMessage[0], str(
                        AtranscriptMessage[1]), str(BtranscriptMessage[1]), delteCode[0], delteCode[1]]
                )+"\n"
                out.append(tmp)

            else:
                pass
    with open(nampeSpace.out, 'w') as File:
        for item in out:
            File.write(item)

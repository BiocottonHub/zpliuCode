'''
Descripttion: 
version: 
Author: zpliu
Date: 2020-12-28 19:44:40
LastEditors: zpliu
LastEditTime: 2020-12-29 20:56:59
@param: 
'''
from getorfbylength import getisoformORF
import sys
import re
from readgtf import getTranscriptInfo
if __name__ == "__main__":
    gtfFile = sys.argv[1]
    genomeFile = sys.argv[2]
    EMBOSSFile = sys.argv[3]
    outFile = sys.argv[4]
    print('reading Isoform Annotion Message...')
    IsoformAnnotion = getTranscriptInfo(gtfFile)
    print('reading Isoform ORF Message...')
    ORFMessage = getisoformORF(gtfFile, EMBOSSFile)
    print('operate Isoform object...')
    with open(outFile, 'w') as File:
        for key in ORFMessage.keys():
            start, end = ORFMessage[key].getAbsoulteORF()
            sequence = IsoformAnnotion[key].getCDSsequence(
                start, end, genomeFile)
            File.write(sequence)

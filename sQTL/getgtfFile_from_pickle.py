'''
Descripttion: 
version: 
Author: zpliu
Date: 2021-03-14 17:23:27
LastEditors: zpliu
LastEditTime: 2021-03-14 17:46:41
@param: 
'''
import sys
import pickle
from readgtf import transcript


def getTranscriptMessage(transcriptDict, outFile):
    out = []
    for key, transcriptObject in transcriptDict.items():
        out.append(transcriptObject.get_transcriptGtf())
        out.append(transcriptObject.get_exonGtf())
    with open(outFile, 'w') as File:
        for item in out:
            File.write(item)


if __name__ == "__main__":
    gff_pickle = sys.argv[1]
    outFile = sys.argv[2]
    with open(gff_pickle, 'rb') as File:
        transcriptDict = pickle.load(File)
    getTranscriptMessage(transcriptDict, outFile)

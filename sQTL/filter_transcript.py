'''
Descripttion: 
version: 
Author: zpliu
Date: 2021-03-11 15:33:51
LastEditors: zpliu
LastEditTime: 2021-03-12 08:44:41
@param: 
'''
import pandas as pd
# ? read gffcompare file
from readgtf import readgtfFile
import pickle
import sys


def getFPKM(AllgeneFPKMFile):
    AllgeneFPKM = pd.read_pickle(AllgeneFPKMFile)
    return AllgeneFPKM


def filter_transcript(gtfFile: str, AllgeneFPKMFile: str, AlltranscriptFile: str):
    # ?load data from pickle
    AllgeneFPKM = getFPKM(AllgeneFPKMFile)
    AlltranscriptFPKM = getFPKM(AlltranscriptFile)
    # ? get transcript object & without intergenic isoforms
    transcriptDict = readgtfFile(gtfFile)
    deltranscriptList = []
    for transcriptID, transcriptObject in transcriptDict.items():
        # ? itera all transcripts id
        geneId = transcriptObject.get_GeneId()
        # ? get FPKM of a transcriptID by row
        geneFPKMrowContent = next(
            AllgeneFPKM.loc[AllgeneFPKM['Gene ID'] == geneId].iterrows())[1]
        transcriptrowContent = next(
            AlltranscriptFPKM.loc[
                AlltranscriptFPKM['t_name'] == transcriptID].iterrows())[1]
        sampleCount = 0
        totalSampleCount = 0
        # todo FPKM below 10% of gene FPKM
        # ?
        for tFPKM, gFPKM in zip(transcriptrowContent[2:], geneFPKMrowContent[1:]):
            totalSampleCount += 1
            if float(gFPKM) == 0:
                sampleCount += 1
            elif float(tFPKM)/float(gFPKM) >= 0.1:
                sampleCount += 1
        # ? sample count below 5%
        if sampleCount/totalSampleCount < 0.05:
            deltranscriptList.append(transcriptID)
    # ! del transcriptId list
    for i in deltranscriptList:
        del transcriptDict[i]
    return transcriptDict


if __name__ == "__main__":
    gtfFile = sys.argv[1]
    AllgeneFPKMpick = sys.argv[2]
    AlltranscriptFPKMpick = sys.argv[3]
    gtf_pickle = sys.argv[4]
    filterTranscript = filter_transcript(
        gtfFile, AllgeneFPKMpick, AlltranscriptFPKMpick)
    with open(gtf_pickle, 'wb') as File:
        pickle.dump(filterTranscript, File)

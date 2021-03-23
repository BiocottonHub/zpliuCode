'''
Descripttion:
version:
Author: zpliu
Date: 2021-03-12 09:05:21
LastEditors: zpliu
LastEditTime: 2021-03-14 08:48:04
@param:
'''
import sys
import pickle
import pandas as pd
from readgtf import transcript


def filter_novel_junctiwion(gtfPickleFile, junctionReadFile):
    with open(gtfPickleFile, 'rb') as f:
        # ? load the gtf file
        transcriptDict = pickle.load(f)
    pass
    # ? junction read count
    junctionRead = pd.read_pickle(junctionReadFile)
    deltranscriptRead = []
    transcriptCount = 0
    for transcriptID, transcriptObject in transcriptDict.items():
        transcriptCount += 1
        #! log the progress
        print(transcriptCount, ":"+transcriptID)
        for chromsome, start, end in transcriptObject.get_intronCoordinate():
            # ? get junction read count in each sample
            sample = 0
            totalsample = 0
            junctionCountList = next(junctionRead.loc[(junctionRead['chr'] == chromsome) & (
                junctionRead['start'] == start) & (junctionRead['end'] == end)].iterrows())[1]
            totalsample += 1
            for readCount in junctionCountList[3:]:
                if readCount >= 10:
                    sample += 1
            # print(sample)
            try:
                if sample/totalsample < 0.05:
                    deltranscriptRead.append(transcriptID)
                    # ? deal with next transcript
                    break
            except ZeroDivisionError:
                #! those transcripts without intron, onle one exon
                pass
    # ? uniq the transcript id
    for item in list(set(deltranscriptRead)):
        del transcriptDict[item]
    #! return the new transcript object
    return transcriptDict


if __name__ == "__main__":
    gffpickleFile = sys.argv[1]
    junctionReadFile = sys.argv[2]
    outpickleFile = sys.argv[3]
    with open(outpickleFile, 'wb') as File:
      #! pickle the object
        pickle.dump(filter_novel_junctiwion(
            gffpickleFile, junctionReadFile), File)

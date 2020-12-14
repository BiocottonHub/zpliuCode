from readFastaFromBed import readFastaFile
from readgtf import getTranscriptInfo
import sys
'''
@ 基因组序列文件
@ 转录本gtf注释文件
@ 输出cDNA序列文件
'''
if __name__ == "__main__":
    genomesequence = readFastaFile(sys.argv[1])  # 读取基因组序列文件
    transcriptDict = getTranscriptInfo(sys.argv[2])  # 读取转录本注释文件
    with open(sys.argv[3], 'w') as File:
        for Key in transcriptDict.keys():
            File.write(transcriptDict[Key].getcDNAsequence(genomesequence))

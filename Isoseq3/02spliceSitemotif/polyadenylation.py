'''
Descripttion: 
version: 
Author: zpliu
Date: 2020-12-30 09:38:49
LastEditors: zpliu
LastEditTime: 2020-12-30 10:15:35
@param: 
'''
import sys
sys.path.insert(0, '/public/home/zpliu/github/zpliuCode/script/genestruct')
while 1:
    from readgtf import getTranscriptInfo
    from readFastaFromBed import getFastaBypybedtools
    break

if __name__ == "__main__":
    gtfFile = sys.argv[1]
    print('reading gtf file')
    transcriptMessage = getTranscriptInfo(gtfFile)
    print('get sequence of polyadenylation...')
    genomeFile = sys.argv[2]

    with open(sys.argv[3], 'w') as File:
        BedtoolsStr = ''
        for transcriptObject in transcriptMessage.values():
            stand = transcriptObject.stand
            exonList = sorted(transcriptObject.getExonCoordinate())
            transcriptId = transcriptObject.transcriptName
            chromsomes = transcriptObject.chromosome
            if stand == "+":
                BedtoolsStr += "\t".join(
                    [chromsomes, str(exonList[-1]-50), str(exonList[-1]), transcriptId, '1', stand])+"\n"
            else:
                BedtoolsStr += "\t".join(
                    [chromsomes, str(exonList[0]-1), str(exonList[0]+49), transcriptId, '1', stand])+"\n"
        File.write(getFastaBypybedtools(BedtoolsStr, genomeFile))

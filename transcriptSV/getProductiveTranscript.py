'''
Descripttion: 
version: 
Author: zpliu
Date: 2020-12-22 19:50:47
LastEditors: zpliu
LastEditTime: 2020-12-22 21:52:27
@param: 
'''
import sys
# this is the var that decide which  column be RNA-seq
DeteeminColumn = 5
RNAseqColumn = 4


def getProdectiveTranscript(transcriptExpression: str):
    '''
    @Descripttion: 
    @param: 
    @return: 
    '''
    out = {}
    with open(transcriptExpression, 'r') as File:
        for line in File:
            line = line.strip("\n").split("\t")
            out[line[0]] = out.get(line[0], [line[1], float(
                line[DeteeminColumn]), float(line[RNAseqColumn])])
            if float(line[DeteeminColumn]) > out[line[0]][1]:
                out[line[0]] = [line[1], float(
                    line[DeteeminColumn]), float(line[RNAseqColumn])]
            elif float(line[DeteeminColumn]) == out[line[0]][1] and float(line[RNAseqColumn]) >= out[line[0]][2]:
                out[line[0]] = [line[1], float(
                    line[DeteeminColumn]), float(line[RNAseqColumn])]
            else:
                pass
    return out


if __name__ == "__main__":
    print(getProdectiveTranscript(sys.argv[1])['Ghir_A01G002010'])

'''
Descripttion: compare the struction of two transcripts,with insert or delet
version: 1.0
Author: zpliu
Date: 2020-12-22 18:57:03
LastEditors: zpliu
LastEditTime: 2020-12-22 20:43:57
@param: 
'''
import os
import re
import sys
from fun.readFasta import readFastaFile


def runMuscle(Atranscript: str, Btranscript: str):
    '''
    @Descripttion: 
    @param: 
    @return: ('0,1,2','2,1,2')
    >Atranscript
    ATACGAAAAAAAAAAGGGGGGGGGGGGGGGGGGGGGCCCCCCCCCCTTTTTTTTTTT
    >Btranscript
    ATACGAAAAAAAAAAGGGGGGGGGGGGGGGGGGGCCCCCCCCCCTTTTTTTTTTT
    #except
    (0,2)
    '''
    Deletpattern = r'[ATCGatcgNn]+(-*)[ATCGatcgNn]+'
    infasta = Atranscript+"\n"+Btranscript+"\n"
    out = os.popen(
        "printf \"%s\"|~/software/muscle3.8.31_i86linux64   2>/dev/null" % infasta).read()  # get a generator
    Aaligint, Baligint = out.split("\n>")
    # 76,45 mean there is 75bp、45bp delet compare with another transctripts
    AdeletCount = ",".join([str(len(i)) for i in re.findall(
        Deletpattern, Aaligint.replace("\n", ""))])  # get the delte count or 0
    BdeletCount = ",".join([str(len(i)) for i in re.findall(
        Deletpattern, Baligint.replace("\n", ""))])
    return (AdeletCount, BdeletCount)


def readFastaFile2(fastaFile: str):
    '''
    @Descripttion:  And file that contained two fasta sequence 
    @param: fastaFile @str
    @return: 
    '''

    out = readFastaFile(fastaFile)
    tmp = ()
    for key, value in out.items():
        tmp += (">"+key+"\n"+value+"\n",)
    return tmp


if __name__ == "__main__":
    A, B = readFastaFile2(sys.argv[1])
    print(runMuscle(A, B))

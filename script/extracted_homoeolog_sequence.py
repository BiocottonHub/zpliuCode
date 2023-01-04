'''
Descripttion: 
version: 
Author: zpliu
Date: 2021-08-02 22:35:21
LastEditors: zpliu
LastEditTime: 2021-08-02 22:53:55
@param: 
'''

from Bio import SeqIO
import pandas as pd
import sys
import re

def readFastaFile(fastaSequenceFile):
    '''read fasta File 
    args: fasta file @str
    return: fasta dict@dict
    '''
    sequenceDict = {}
    for seq_record in SeqIO.parse(open(fastaSequenceFile, mode='r'), 'fasta'):
        geneId=re.split("\.[0-9]*",seq_record.id)[0]
        #geneId = seq_record.id.split("\.")[0]
        if geneId not in sequenceDict:
            sequenceDict[geneId] = seq_record.seq
        elif len(seq_record.seq) > len(sequenceDict[geneId]):
            #! the longest transcripts
            sequenceDict[geneId] = seq_record.seq
        else:
            pass
    return sequenceDict
if __name__=="__main__":
    geneId=pd.read_csv(sys.argv[1],header=None,index_col=None,sep="\t")
    #!  read All gene sequence File
    sequenceDict=readFastaFile(sys.argv[2])
    #! save to file
    with open(sys.argv[3],'w') as File:
        for value in geneId.values:
            geneA,geneB=value

            File.write(
                ">"+geneA+"\n"+"".join(sequenceDict[geneA])+"\n"
            )
            File.write(
                ">"+geneB+"\n"+"".join(sequenceDict[geneB])+"\n"
            )

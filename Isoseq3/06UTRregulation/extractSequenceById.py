'''
Descripttion: 
version: 
Author: zpliu
Date: 2021-01-02 15:14:09
LastEditors: zpliu
LastEditTime: 2021-01-02 15:27:22
@param: fasta file 
@param: gene Id 
@param: out file 
'''
import sys
sys.path.insert(0, '/public/home/zpliu/github/zpliuCode/script/genestruct')

while 1:
    #from readgtf import getTranscriptInfo
    from readFastaFromBed import readFastaFile
    break


if __name__ == "__main__":
    fastaDict = readFastaFile(sys.argv[1])
    with open(sys.argv[3], 'w') as File:
        with open(sys.argv[2], 'r') as File1:
            for line in File1:
                line = line.strip("\n").split("\t")
                try:
                    File.write(">"+line[0]+"\n"+fastaDict[line[0]]+"\n")
                except KeyError:
                    print("WArning:"+line[0]+"\tnot in the fasta sequence")

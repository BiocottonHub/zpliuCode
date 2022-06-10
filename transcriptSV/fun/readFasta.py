'''
Descripttion: 
version: 
Author: zpliu
Date: 2020-12-22 19:41:30
LastEditors: zpliu
LastEditTime: 2020-12-23 21:57:42
@param: 
'''
import re
import sys


def readFastaFile(fastaFile: str):
    '''
    @Descripttion:  And file that contained two fasta sequence 
    @param: fastaFile @str
    @return: 
    '''
    out = {}
    tmp = ''
    with open(fastaFile, 'r') as File:
        for line in File:
            if re.match("^>", line):
                tmp = line.strip("\n").strip(">").split("\s+")[0]
            else:
                out[tmp] = out.get(tmp, '')+line.strip("\n")
    return out


if __name__ == "__main__":
    print(readFastaFile(sys.argv[1]))

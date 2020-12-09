'''
Descripttion: 
version: 
Author: zpliu
Date: 2020-12-07 09:55:47
LastEditors: zpliu
LastEditTime: 2020-12-08 20:00:04
@param: 
'''
import sys
import gzip
import re


def extract_read_count(inputfile: tuple):
    '''
    @Descripttion: extract read Count from raw file
    @param: CpG|CHG|CHH context file @str
    @return: out file @ dict
    '''
    out = {}
    with open(inputfile, 'r') as File:
        line = File.readline().strip("\n").split("\t")
        if len(line) == 5:
            # raw methylation file
            if line[1] == "+":
                out[line[0]+"*"+line[1]] = [1, 0]
            else:
                out[line[0]+"*"+line[1]] = [0, 1]
            for line in File:
                line = line.strip("\n").split("\t")
                if line[1] == "+":
                    out[line[2]+"*"+line[3]
                        ] = out.get(line[2]+"*"+line[3], [0, 0])
                    out[line[2]+"*"+line[3]][0] += 1
                else:
                    out[line[2]+"*"+line[3]
                        ] = out.get(line[2]+"*"+line[3], [0, 0])
                    out[line[2]+"*"+line[3]][1] += 1
        else:
            # split file
            out[line[0]+"*"+line[1]] = out.get(line[0]+"*"+line[1], [0, 0])
            out[line[0]+"*"+line[1]][0] += int(line[2])
            out[line[0]+"*"+line[1]][1] += int(line[3])
            for line in File:
                line = line.strip("\n").split("\t")
                out[line[0]+"*"+line[1]] = out.get(line[0]+"*"+line[1], [0, 0])
                out[line[0]+"*"+line[1]][0] += int(line[2])
                out[line[0]+"*"+line[1]][1] += int(line[3])
    return out


if __name__ == "__main__":
    out = extract_read_count(sys.argv[1])
    with open(sys.argv[2], 'w') as File:
        for key, value in out.items():
            outStr = "\t".join(key.split("*"))+"\t" + \
                str(value[0])+"\t"+str(value[1])+"\n"
            File.write(outStr)

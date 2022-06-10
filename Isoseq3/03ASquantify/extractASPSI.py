'''
Descripttion: 
version: 
Author: zpliu
Date: 2021-01-04 21:53:23
LastEditors: zpliu
LastEditTime: 2021-01-05 15:06:52
@param: 
'''
import sys


def caculatePSI(PSIfile):
    out = {}
    with open(PSIfile, 'r') as File:
        for line in File:
            line = line.strip("\n").split("\t")
            if(int(line[2]) == 0):
                out[line[1]] = '0'+"\t"+str((int(line[2])+int(line[3])))
            elif(int(line[3]) == 0):
                out[line[1]] = '1'+"\t"+str((int(line[2])+int(line[3])))
            else:
                out[line[1]] = str(
                    int(line[2])/(int(line[2])+int(line[3])))+"\t"+str((int(line[2])+int(line[3])))
    return out


if __name__ == "__main__":
    out = []
    ASPSI = caculatePSI(sys.argv[1])
    with open(sys.argv[2], 'r') as File1:
        for line in File1:
            line = line.strip("\n").split("\t")
            out.append("\t".join(line)+"\t" +
                       str(ASPSI[line[0]])"\n")
    with open(sys.argv[3], 'w') as File:
        for item in out:
            File.write(item)

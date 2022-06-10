'''
Descripttion: 
version: 
Author: zpliu
Date: 2021-02-07 21:16:27
LastEditors: zpliu
LastEditTime: 2021-02-09 15:48:41
@param: 
'''
import sys
if __name__ == "__main__":
    out = {}
    with open(sys.argv[1], 'r') as File:
        for line in File:
            line = line.strip("\n").split("\t")
            out[line[0]] = out.get(line[0], [0, 0])
            out[line[0]][0] += int(line[1])
            out[line[0]][1] += int(line[2])
    with open(sys.argv[1], 'w') as File:
        for key, values in out.items():
            try:
                File.write(key+"\t"+str(values[0]/values[1])+"\n")
            except ZeroDivisionError:
                File.write(key+"\t"+str(0)+"\n")

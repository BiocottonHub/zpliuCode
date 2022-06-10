'''
Descripttion: 
version: 
Author: zpliu
Date: 2021-01-06 09:31:58
LastEditors: zpliu
LastEditTime: 2021-01-06 09:43:13
@param: 
'''
import sys
if __name__ == "__main__":
    out = {}
    with open(sys.argv[1], 'r') as File:
        for line in File:
            line = line.strip("\n").split("\t")
            out[line[0]] = out.get(
                line[0], [line[1], int(line[2]), int(line[3])])
            if abs(int(line[2])-int(line[3])) < abs(out[line[0]][1]-out[line[0]][2]):
                out[line[0]] = [line[1], int(line[2]), int(line[3])]
            else:
                pass
    out2 = {}
    for key, value in out.items():
        if value[0] not in out2:
            out2[value[0]] = out2.get(value[0], [key, value[1], value[2]])
        else:
            if abs(value[1]-value[2]) < abs(out2[value[0]][1]-out2[value[0]][2]):
                out2[value[0]] = [key, value[1], value[2]]
    with open(sys.argv[2], 'w') as File:
        for key, value in out2.items():
            File.write(
                value[0]+"\t"+key+"\t"+str(value[1])+"\t"+str(value[2])+"\n"
            )

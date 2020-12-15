'''
Descripttion: 
version: 
Author: zpliu
Date: 2020-11-18 15:51:47
LastEditors: zpliu
LastEditTime: 2020-11-18 16:30:57
@param: 
'''
import re
import sys
import argparse


def readPfamFile(PfamFile):
    # 获取每个isoform的所有保守结构域
    out = {}
    with open(PfamFile, 'r') as File:
        for line in File:
            if not re.match(r'^[#\s+]', line):  # 去除注释行和空白行
                line = re.split(r'\s+', line.strip("\n"))
                out[line[0]] = out.get(line[0], {
                    line[5]: []
                })
                if line[5] not in out[line[0]]:
                    out[line[0]][line[5]] = [str(line[3])+"-"+str(line[4])]
                else:
                    out[line[0]][line[5]].append(str(line[3])+"-"+str(line[4]))
    for key in out:
        for key2 in out[key]:
            out[key][key2].sort()  # 给所有的保守结构域排序数组
    return out


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="")
    parser.add_argument("-APfam", help="reference ORF")
    parser.add_argument("-BPfam", help="reference ORF")
    parser.add_argument("-conserve", help="reference ORF")
    parser.add_argument("-o", help="out put File")
    args = parser.parse_args()
    Apfam = readPfamFile(args.APfam)
    Bpfam = readPfamFile(args.BPfam)

    out = []
    with open(args.conserve, 'r') as File:
        for line in File:
            line = line.strip("\n").split("\t")
            if line[3] != line[7]:
                continue
            elif line[3] == "None":
                continue
            else:
                if Apfam[line[0]] == Bpfam[line[4]]:
                    out.append("\t".join(line)+"\n")
    with open(args.o, 'w') as File:
        for line in out:
            File.write(line)

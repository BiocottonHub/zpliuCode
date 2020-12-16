'''
Descripttion: 
version: 
Author: zpliu
Date: 2020-12-01 20:46:02
LastEditors: zpliu
LastEditTime: 2020-12-01 23:06:11
@param: -a 
@param: -b
@param: -o
-a : scale chromsome file
Chr01	0	113035596	Ghir_A01	0	117710661
Chr02	0	99090824	Ghir_A02	0	108049532
Chr03	0	135709677	Ghir_A03	0	113014280
Chr04	0	98600468	Ghir_A04	85114396	0
Chr05	0	97880472	Ghir_A05	0	109365994
'''
import argparse
import re


def deteminator(chromsomes, length1, length2, inverted):
    '''
    @Descripttion: 
    @param: 
    @return: 
    @ 这里不存在闭包问题，因为最外层函数每执行一次才返回一个子函数
    '''
    scaled = length2/length1

    def function1(start, end):
        return chromsomes+"\t"+str(int(start*scaled))+"\t"+str(int(end*scaled))

    def function2(start, end,):
        end2 = ((length1-start)/length1)*length2
        start2 = ((length1-end)/length1)*length2
        return chromsomes+"\t"+str(int(start2))+"\t"+str(int(end2))

    if inverted:
        return function2
    else:
        return function1


def readScaleChromsome(chromsomeFile):
    out = {}
    with open(chromsomeFile, 'r') as File:
        for line in File:
            line = line.strip("\n").split("\t")
            if int(line[-1]) > int(line[-2]):
                out[line[0]] = deteminator(
                    line[3], int(line[2]), int(line[5]), 0)
            else:
                out[line[0]] = deteminator(
                    line[3], int(line[2]), int(line[4]), 1)
    return out


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-a", help="scaled chromsomed bed")
    parser.add_argument("-b", help="feature of A genome bed ")
    parser.add_argument("-o", help="scale A genome to B genome bed")
    args = parser.parse_args()
    scaletor = readScaleChromsome(args.a)
    out = []
    with open(args.b, 'r') as File:
        for line in File:
            line = line.strip("\n").split("\t")
            bed = scaletor[line[0]](int(line[1]), int(line[2]))
            if len(line)>3:
                out.append(bed+"\t"+"\t".join(line[3::])+"\n")
            else:
                out.append(bed+"\n")
    with open(args.o, 'w') as File:
        for item in out:
            File.write(item)

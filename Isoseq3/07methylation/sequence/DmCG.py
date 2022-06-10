'''
Descripttion: 
version: 
Author: zpliu
Date: 2021-02-19 10:26:32
LastEditors: zpliu
LastEditTime: 2021-02-19 14:39:39
@param: 
'''
import sys
import os
from scipy import stats


def getSite(locationStr):
    return locationStr.split("-")[1]


def getMethylationRatio(Location, MethylationFile):
    MethlationLines = os.popen(
        'samtools view -O SAM %s %s' % (MethylationFile, Location)).read()
    if MethlationLines == "":
        return 'None'  # none CG
    else:
        for item in MethlationLines.split("\n")[:-1]:
            MethlationLocal = getSite(item.split("\t")[0])
            if MethlationLocal == getSite(Location):  # same C Site
                read1 = item.split("\t")[0].split("-")[2]
                read2 = item.split("\t")[0].split("-")[3]
                if eval(read1+"+"+read2) >= 3:
                    return str(eval(read1+"/("+read1+"+"+read2+")"))
                else:
                    return '0'  # less than 3 read
            else:
                pass
    return 'None'  # none same C site


def getPvalue(rep1, rep2, rep3, rep4):
    def f(x): return 0 if x == "None" else float(x)
    mean1 = (f(rep1)+f(rep2))/2
    mean2 = (f(rep3)+f(rep4))/2
    if abs(mean1-mean2) < 0.5:
        return 'None'
    else:
        F, p = stats.f_oneway((f(rep1), f(rep2)), (f(rep3), f(rep4)))
        return str(p)


def getDmCG(LocationFile, MethylationFileR1, MethylationFileR2, MethylationFileR3, MethylationFileR4, outFile):
    out = []
    with open(LocationFile, 'r') as File:
        for line in File:
            line = line.strip("\n").split("\t")
            location1 = line[0]
            location2 = line[2]
            rep1 = getMethylationRatio(location1, MethylationFileR1)
            rep2 = getMethylationRatio(location1, MethylationFileR2)
            rep3 = getMethylationRatio(location2, MethylationFileR3)
            rep4 = getMethylationRatio(location2, MethylationFileR4)
            if rep1 == "None" and rep2 == 'None' and rep3 == 'None' and rep4 == "None":
                pass
            else:
                pvalue = getPvalue(rep1, rep2, rep3, rep4)
                tmp = "\t".join(line)+"\t" + \
                    "\t".join((rep1, rep2, rep3, rep4, pvalue))+"\n"
                print(tmp)
                out.append(tmp)
    with open(outFile, 'w') as File:
        for item in out:
            File.write(item)


if __name__ == "__main__":
    LocationFile = sys.argv[1]
    outFile = sys.argv[2]
    MethylationFileR1 = sys.argv[3]
    MethylationFileR2 = sys.argv[4]
    MethylationFileR3 = sys.argv[5]
    MethylationFileR4 = sys.argv[6]
    getDmCG(LocationFile, MethylationFileR1, MethylationFileR2,
            MethylationFileR3, MethylationFileR4, outFile)

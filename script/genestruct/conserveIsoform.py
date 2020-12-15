
import re
import sys
import argparse


def readBlastFile(blastFile):
    '''
    根据转录本CDS blastp的结果；
    判断两个转录本间是否存在相似度大于90%的片段

    '''
    out = {}
    with open(blastFile, 'r') as File:
        for line in File:
            line = line.split("\t")
            identity = float(line[-4])
            if identity >= 90:
                if line[0] not in out:
                    out[line[0]] = {}
                    out[line[0]][line[1]] = line[-4]
                else:
                    if line[1] not in out[line[0]]:
                        out[line[0]][line[1]] = line[-4]
                    else:
                        pass
    return out


def readisoformMessage(isoformFile):
    '''
    isoform的信息还包括转录本的长度信息
    '''
    out = {}
    with open(isoformFile, 'r') as File:
        for line in File:
            line = line.strip("\n").split("\t")
            if line[0] not in out:
                out[line[0]] = [{
                    'isoformName': line[1],
                    'isoformLength':line[2]+"\t"+line[3]
                }]
            else:
                out[line[0]].append({
                    'isoformName': line[1],
                    'isoformLength': line[2]+"\t"+line[3]
                })
    return out


def readhomologGene(homologGeneFile, column1, column2):
    out = {}
    with open(homologGeneFile, 'r') as File:
        for line in File:
            line = line.strip("\n").split("\t")
            out[line[column1-1]] = line[column2-1]
            #out[line[column2-1]] = line[column1-1]
    return out  # 两个同源基因组成的关联字典


def readPfamFile(PfamFile, genome):
    # 获取每个isoform的所有保守结构域
    out = {}
    with open(PfamFile, 'r') as File:
        for line in File:
            if not re.match(r'^[#\s+]', line):  # 去除注释行和空白行
                line = re.split(r'\s+', line.strip("\n"))
                line[0] = genome+line[0]
                if line[0] not in out:
                    out[line[0]] = [line[5]]
                else:
                    out[line[0]].append(line[5])
    for key in out:
        out[key].sort()  # 给所有的保守结构域排序数组
    return out


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="")
    parser.add_argument("-Aisoform", help="reference gtf file")
    parser.add_argument("-Bisoform", help="PacBio gtf file")
    parser.add_argument("-APfam", help="reference ORF")
    parser.add_argument("-BPfam", help="reference ORF")
    parser.add_argument("-blast", help="reference ORF")
    parser.add_argument("-homolog", help="PacBio ORF")
    parser.add_argument("-o", help="out put File")
    args = parser.parse_args()
    homologDict = readhomologGene(args.homolog, 1, 6)  # 同源基因信息
    Aisoform = readisoformMessage(args.Aisoform)
    Bisoform = readisoformMessage(args.Bisoform)
    APfam = readPfamFile(args.APfam, 'At-')
    BPfam = readPfamFile(args.BPfam, 'Dt-')
    blastDict = readBlastFile(args.blast)
    # 遍历所有的同源基因信息
    out = []
    for key in homologDict:
        # 主要是看基因有没有表达PacBio转录本
        if key in Aisoform and homologDict[key] in Bisoform:
            isoformsA = Aisoform[key]  # 获取基因的所有PacBio isofrom信息
            isoformsB = Bisoform[homologDict[key]]
            for item1 in isoformsA:
                for item2 in isoformsB:
                    # blast相似度大于90%
                    try:
                        identity = blastDict[item1['isoformName']
                                             ][item2['isoformName']]
                    except KeyError:
                        identity = '0'

                    if item1['isoformName'] in APfam and item2['isoformName'] in BPfam:
                        out.append(item1['isoformName']+"\t" +
                                   item1['isoformLength']+"\t"+",".join(APfam[item1['isoformName']]) +
                                   "\t"+item2['isoformName']+"\t"+item2['isoformLength']+"\t" +
                                   ",".join(BPfam[item2['isoformName']]) +
                                   "\t"+identity+"\n")
                    elif item1['isoformName'] in APfam and item2['isoformName'] not in BPfam:
                        out.append(item1['isoformName']+"\t" +
                                   item1['isoformLength']+"\t"+",".join(APfam[item1['isoformName']]) +
                                   "\t"+item2['isoformName']+"\t"+item2['isoformLength']+"\t" +
                                   "None\t"+identity+"\n")
                    elif item1['isoformName'] not in APfam and item2['isoformName'] in BPfam:
                        out.append(item1['isoformName']+"\t" +
                                   item1['isoformLength']+"\tNone" +
                                   "\t"+item2['isoformName']+"\t"+item2['isoformLength']+"\t" +
                                   ",".join(BPfam[item2['isoformName']]) +
                                   "\t"+identity+"\n")
                    else:
                        out.append(item1['isoformName']+"\t" +
                                   item1['isoformLength']+"\tNone" +
                                   "\t"+item2['isoformName']+"\t"+item2['isoformLength']+"\t" +
                                   "None\t"+identity+"\n")

    with open(args.o, 'w') as File:
        for item in out:
            File.write(item)

from readgtf import getGeneInfo
import argparse
# 创建参数对象
parser = argparse.ArgumentParser(description="")
parser.add_argument("-refgtf", help="At blast out file")
parser.add_argument("-isogtf", help="Dt blast out file")
parser.add_argument("-o", help="D5 blast out file")
# 解析参数对象
args = parser.parse_args()
refGene = getGeneInfo(args.refgtf)
PacBiogene = getGeneInfo(args.isogtf)
out = []
for gene in PacBiogene.keys():
    allsplitesite = refGene[gene].getAllspliceSite()
    for item in PacBiogene[gene].getAlltranscriptsObject():
        intronCoordinate = item.getIntronCoordinate()
        flag = True
        a = 0  # 与参考基因组一直的剪切位点数
        b = 0  # 与参考基因组不一直的剪切位点数
        for intron in intronCoordinate:  # 将每个剪切位点与参考基因组的剪切位点进行比较
            if intron[0] not in allsplitesite or intron[1] not in allsplitesite:
                flag = False
                b += 1
            else:
                a += 1
                continue
        if flag:
            out.append(item.getTranscriptName()+"\t"+str(item.getcDNAlength())+"\t"+gene +
                       "\tAnnotion\t"+str(a)+"\t"+str(b)+"\n")
        else:
            out.append(item.getTranscriptName()+"\t"+str(item.getcDNAlength())+"\t"+gene +
                       "\tunAnnotion\t"+str(a)+"\t"+str(b)+"\n")
            # print(item.getTranscriptName()+"\t"+gene)
            # print(intronCoordinate)
            # print(allsplitesite)
            # breakpoint


with open(args.o, 'w') as File:
    for item in out:
        File.write(item)

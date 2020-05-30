'''/*
 * @Author: zpliu
 * @Date: 2020-05-29 15:40:14
 * @Last Modified by:   zpliu
 * @Last Modified time: 2020-05-29 15:40:14
 *Usage:
 *python script.py -h
 * -t  迭代次数，生成的文件数
 * -d 输出目录
 * -TAD TAD的最小长度
 * -init 初始bodary 随机位置
 */
'''
import argparse
import random
import os
parser = argparse.ArgumentParser(
    description="intersect conserve  homolog AS event")
parser.add_argument("-chr", help="chrosomoe bed file")
parser.add_argument("-t", help="lterate times")
parser.add_argument("-d", help="out directory")
parser.add_argument("-TAD", help="the minimal length of TAD")
parser.add_argument("-init", help="init boundary random position")

args = parser.parse_args()


def mkdir(path):
    # 去除首位空格
    path = path.strip()
    # 去除尾部 \ 符号
    path = path.rstrip("\\")
    isExists = os.path.exists(path)
    if not isExists:
        os.makedirs(path)
        print(path+' 创建成功')
        return True
    else:
        # 如果目录存在则不创建，并提示目录已存在
        print(path+' 目录已存在')
        return False


def satisfGap(Array1, gap):
    for index in range(1, len(Array1)):
        if Array1[index] - Array1[index-1] < gap:
            return False
        else:
            continue
    return True


def ifoverchromosome(anchro, chromosomelength):
    '''
    判断bodunary位置与染色长度间关系
    '''
    if anchro+int(args.gap) < chromosomelength and anchro-int(args.gap) > 1:
        return 1
    if anchro+int(args.gap) < chromosomelength and anchro-int(args.gap) < 1:
        return 2
    if anchro+int(args.gap) > chromosomelength and anchro-int(args.gap) > 1:
        return 3


def getAnchroNum(anchro):
    '''
    对指定范围的染色体获取随机数组，
    再从随机数组中获取随机数
    '''
    tmp = []
    for item in anchro:
        tmp.append(random.randint(item[0], item[1]))
    return random.choice(tmp)


def intersectAnchro(anchro1, anchro2):
    '''
    每添加一个boundary就更新，锚点区域
    '''
    tmp = []
    for item1 in anchro1:
        for item2 in anchro2:
            if item1[0] > item2[1] or item1[1] < item2[0]:
                continue
            else:
                tmp2 = item1+item2
                tmp2.sort()
                tmp.append([tmp2[1], tmp2[2]])
    if(len(tmp) == 0):
        print("锚点出错!!!!!!!!!!!!!!")  # 取交集时出现了问题
    else:
        return tmp


mkdir(args.d)
out = []
with open(args.chr, 'r') as File:
    for line in File.readlines():
        line = line.strip("\n").split("\t")
        chromosomeName = line[0]
        end = int(line[1])
        count = int(line[2])+1
        mkdir(args.d+"/"+chromosomeName)
        i = 1
        while(i <= int(args.t)):
            with open(args.d+"/"+chromosomeName+"/"+str(i), 'w') as outFile:
                while(len(out) < count):
                    if len(out) == 0:  # 获取初始boundary的随机坐标
                        out.append(random.randint(1, int(args.init)))
                    elif len(out) == 1:  # 其实这步可以和else合并
                        if ifoverchromosome(out[0], end) == 1:
                            out.append(getAnchroNum(
                                [[1, out[0]-int(args.gap)], [out[0]+int(args.gap), end]]))
                        elif ifoverchromosome(out[0], end) == 2:
                            out.append(getAnchroNum(
                                [[out[0]+int(args.gap), end]]))
                        elif ifoverchromosome(out[0], end) == 3:
                            out.append(getAnchroNum(
                                [[1, out[0]-int(args.gap)]]))
                    else:  # 对多个boundary进行计算
                        out.sort()
                        anchroindex = [[1, end]]
                        for item in out:
                            if ifoverchromosome(item, end) == 1:
                                anchroindex = intersectAnchro(
                                    anchroindex, [[1, item-int(args.gap)], [item+int(args.gap), end]])
                            elif ifoverchromosome(item, end) == 2:
                                anchroindex = intersectAnchro(
                                    anchroindex, [[item+int(args.gap), end]])
                            elif ifoverchromosome(item, end) == 3:
                                anchroindex = intersectAnchro(
                                    anchroindex, [[1, item-int(args.gap)]])
                        out.append(getAnchroNum(
                            anchroindex))
                outFile.write("\n".join([str(i) for i in out]))
                out = []
                print(chromosomeName+"\t的第"+str(i)+"\t文件完成")
                i += 1

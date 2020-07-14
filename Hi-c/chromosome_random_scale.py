# encoding:utf-8
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
import numpy as np
import os
import sys
import multiprocessing
parser = argparse.ArgumentParser(
    description="intersect conserve  homolog AS event")
parser.add_argument("-chr", help="chrosomoe bed file")
parser.add_argument("-t", help="lterate times")
parser.add_argument("-d", help="out directory")
parser.add_argument("-p", help="number of threads")
parser.add_argument("-min", help="the minimal length of TAD")
parser.add_argument("-max", help="the minimal length of TAD")

args = parser.parse_args()

'''
随机得到TAD的长度
'''


def getrandomTADlength():
    return random.randint(int(args.min), int(args.max))


'''
截取TAD后，剩余染色体坐标
TADarray [start, end]
chrosomeArray:[[start1,end1],[start2,end2]]
'''
def remaindChrosome(TADarray, chrosomeArray):
    for i in chrosomeArray:
        if i[0] > TADarray[1] and i[1] < TADarray[0]:  # 与染色体片段没有交集
            continue
        else:  # 找到TAD所在的染色体区域
            break
    chrosomeArray.remove(i)  # 去除TAD所在区域
    if TADarray[0]-i[0] >= int(args.min):
        chrosomeArray.append([i[0], TADarray[0]-1])
    if i[1]-TADarray[1] >= int(args.min):
        chrosomeArray.append([TADarray[1]+1, i[1]])
    return chrosomeArray

'''
根据文件夹的名字创建目录
'''
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

'''
重指定的染色体范围随机出一定长度的TAD, 有可能随机不出结果陷入死循环
'''
def randomTAD(TADlength, chrosomeArray):
    tmp = []
    for i in chrosomeArray:
        if i[1]-i[0] < TADlength:
            continue
        else:
            start = np.arange(i[0], i[1]-TADlength, 1)
            tmp += [i[0], i[1]-TADlength]
    #
    tmp2 = []
    for i in range(0, len(tmp), 2):
        # 获取每个片段的随机起始位置
        tmp2.append(random.sample(range(tmp[i], tmp[i+1]), 1)[0])
    try:
        start = random.choice(tmp2)  # 从随机起始位置挑一个
    except:
        print("找不到满足的TAD了,重新获取一个随机TAD长度")
        return []
    return [start, start+TADlength-1]

'''
###使用递归来随机TAD的长度信息
# total 染色体内需要找到的TAD数目
# count 记录成功迭代的数目
# chrosomeArray 记录除去TAD后剩余的长度
'''
def RecursicvegetTAD(total, count, chrosomeArray):
    if count == total:
        return []
    else:
        # 获取本次随机TAD长度
        TADlength = getrandomTADlength()
        # 获取TAD坐标数组
        i = 100
        while i >= 0:  # 循环100次，防止得不到TAD
            TADArray = randomTAD(TADlength, chrosomeArray)
            if(len(TADArray)) != 0:
                break
            i = i-1
        if i < 0:  # 最后一次循环都没搞出有效的TAD，这轮迭代作废
            return [-1]
        # print("找到第"+str(count+1)+"随机TAD,进入下一层递归")

        chrosomeArray = remaindChrosome(TADArray, chrosomeArray)
        # print(chrosomeArray)
        count += 1
        return TADArray+RecursicvegetTAD(total, count, chrosomeArray)

'''
每个进程执行的函数
传入需要执行的染色体数据，和进程ID号
每条染色体体进行1000轮迭代，完成一轮迭代就把一个结果写入文件
##如果某一轮迭代进入死循环，则重新本轮迭代
'''
def processFile(lines, i):
    print(">>>进程:"+str(i)+"启动将处理:" +
          "\t".join([i.split("\t")[0] for i in lines]))
    for line in lines:
        line = line.strip("\n").split("\t")
        mkdir(args.d+"/"+line[0])
        if(int(line[2]) % 2 == 0):  # boundary为奇数的情况
            total = int(line[2])/2
        else:
            total = int(line[2])/2+1
        for i in range(0, int(args.t)):
            out = RecursicvegetTAD(total, 0, [[0, int(line[-2])]])
            while out[-1] == -1:  # 这轮迭代陷入死循环，作废重新迭代
                print("第："+str(i)+"轮迭代陷入死循环，重新本轮迭代")
                out = RecursicvegetTAD(total, 0, [[0, int(line[-2])]])
            print("染色体:"+line[0]+"完成第"+str(i+1)+"轮迭代")
            with open(args.d+"/"+line[0]+"/"+str(i+1), 'w') as File2:
                for item in out:
                    File2.write(line[0]+"\t"+str(item-100000) +
                                "\t"+str(item+100000)+"\n")


with open(args.chr, 'r') as File:
    lines = File.readlines()
# processFile(lines)


p = multiprocessing.Pool(20)
average = int(len(lines)/int(args.p))
for i in range(0, int(args.p)):
    if i == int(args.p)-1:
        start = i*average
        end = len(lines)
    else:
        start = i*average
        end = (i+1)*average
    p.apply_async(processFile, (lines[start:end], i),
                  callback=None)  # 传参数时，逗号不能够少
p.close()
p.join()

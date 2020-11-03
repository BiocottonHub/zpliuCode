#!/public/home/software/opt/bio/software/Python/2.7.15/bin/python
# coding=utf-8
'''/*
 * @Author: zpliu 
 * @Date: 2020-10-20 15:54:07 
 * @Last Modified by:   zpliu 
 * @Last Modified time: 2020-10-20 15:54:07 
 * @Usage: module load peer/1.0 & python2 peer_interface.py -h 
 */
 '''

import peer
import scipy as SP
import numpy as np
import os
import argparse


def mkdir(path):
    path = path.strip()
    # 去除尾部 \ 符号
    path = path.rstrip("\\")
    # 判断路径是否存在
    # 存在     True
    # 不存在   False
    isExists = os.path.exists(path)
    # 判断结果
    if not isExists:
        os.makedirs(path)
        print(path+' 创建成功')
        return True
    else:
        # 如果目录存在则不创建，并提示目录已存在
        print(path+' 目录已存在')
        return False


def peerModel(expressFile, factorNum):
    model = peer.PEER()
    model.setPhenoMean(expressFile)
    model.setNk(factorNum)
    model.update()
    return model


def saveResiduals(dirName, model):
    residuals = model.getResiduals()
    np.savetxt(dirName+"/residuals.csv", residuals, fmt='%f', delimiter=',')


def saveinferredFactors(dirName, model):
    factors = model.getX()
    np.savetxt(dirName+"/X.csv", factors, fmt='%f', delimiter=',')


def saveweightsFactors(dirName, model):
    weights = model.getW()
    np.savetxt(dirName+"/W.csv", weights, fmt='%f', delimiter=',')


def saveinverseVariance(dirName, model):
    precision = model.getAlpha()
    np.savetxt(dirName+"/Alpha.csv", precision, fmt='%f', delimiter=',')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="")
    parser.add_argument("-f", help="expression file")
    parser.add_argument("-n", help="number of factors")
    parser.add_argument("-has_rownames", help="with header ")
    parser.add_argument("-has_header", help="with rowNames")
    parser.add_argument("-o", help="out put dir")
    args = parser.parse_args()
    '''
    创建输出目标
    '''
    if not mkdir(args.o):
        exit()
    else:
        pass
    '''
    读取文件时，判断是否存在行名和列名
    '''
    if args.has_header:
        expr = SP.loadtxt(args.f, delimiter=',', skiprows=1)
    else:
        expr = SP.loadtxt(args.f, delimiter=',')
    if args.has_rownames:
        expr = expr[:, 1:]
    else:
        pass
    # 模型计算
    model = peerModel(expr, int(args.n))
    # 保存文件
    saveResiduals(args.o, model)
    saveinferredFactors(args.o, model)
    saveweightsFactors(args.o, model)
    saveinverseVariance(args.o, model)

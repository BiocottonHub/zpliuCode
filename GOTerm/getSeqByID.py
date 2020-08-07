'''/*
 * @Author: zpliu 
 * @Date: 2020-08-07 09:31:31 
 * @Last Modified by:   zpliu 
 * @Last Modified time: 2020-08-07 09:31:31 
 * 根据基因编号与物种编号获取对应的fasta序列信息
 * 得到的序列信息进行blast
 */
'''
import urllib.error
import urllib.request
import urllib.parse
import json
import time
import sys
import re


def getSeq(geneName, taxonId):
    # 匹配基因起始和终止坐标
    URLpattern = 'a title="Nucleotide FASTA report" href="([^"]*)"'
    getGoItemUrl = "https://www.ncbi.nlm.nih.gov/gene?term=("+geneName + \
        "%5BGene%20Name%5D)%20AND%20"+taxonId+"%5BTaxonomy%20ID%5D"
    response = urllib.request.urlopen(
        getGoItemUrl).read().decode('utf-8')  # 发起get请求
    # print(response.read().decode('utf-8'))
    try:
        getFastaURL = re.search(URLpattern, response).group(1)  # 获取基因的坐标
    except AttributeError:
        print("\033[1;31m 当前基因编号: "+geneName +
              " 物种编号: "+taxonId+"在数据库中未找到\033[0m!")  # 终端显示红色
        return False
    '''
    获取NCNI中物种编号和基因坐标
    '/nuccore/NC_030079.1?report=fasta&amp;from=38409546&amp;to=38415063'
    '''
    # 获取物种基因组版本编号<img usemap="#nbrMap821891_240255695"
    taxonId = re.search(
        r'<img usemap="([^"]*)"', response).group(1).split("_")[1]  # 获取基因组编号
    start = re.search(r'from=([0-9]*)\&', getFastaURL).group(1)  # 获取开始位点
    end = re.search(r'to=([0-9]*)', getFastaURL).group(1)  # 获取结束位点
    getFastaURL = 'https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id='+taxonId + \
        '&db=nuccore&report=fasta&retmode=text&withmarkup=on&tool=portal&log$=seqview&from=' + \
        start+'&to='+end+'&strand=on'
    # 爬取基因fasta序列信息,统一为正链序列，而不是有义链序列
    response = urllib.request.urlopen(getFastaURL)
    fasta = response.read().decode('utf-8')
    return fasta


outSeq = []
'''
爬取GO时，获取对应的列
基因编号: 第7列
物种编号: 物种编号
'''
with open(sys.argv[1], 'r') as File:
    for line in File.readlines():
        line = line.strip("\n").split("\t")
        try:
            taxonId = line[10].split(":")[1]  # 出错就是表头的问题
        except:
            continue
        fasta = getSeq(line[6], taxonId)
        if fasta:
            print("爬取基因序列信息: "+line[6])
            outSeq.append(fasta+"\n",)
        else:
            continue

with open(sys.argv[2], 'w') as File:
    for item in outSeq:
        File.write(item)

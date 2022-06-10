# -*- coding:UTF-8 -*-
#
'''
author:zpliu
date:2019-07-08
email:1944532210@qq.com
descripte: Thie script is use for batch request from 
http://gene-regulation.com/cgi-bin/pub/programs/patch/bin/patch.cgi;and parse the html content into text
usage: python3 patch.py geneSequence.fasta output.txt
'''
import urllib.request
import urllib.parse
import http.cookiejar
import time
import sys
from bs4 import BeautifulSoup
from tqdm import tqdm
import os
# 导入读取基因文件的函数


def fastaread(fastaFile):
    out = {}
    with open(fastaFile, 'r') as File:
        for line in File:
            line = line.strip("\n")
            if line.startswith(">"):
                geneId = line.strip(">")
                out[geneId] = out.get(geneId, '')
            else:
                out[geneId] += line
    return out


# get_url为使用cookie所登陆的网址，该网址必须先登录才可
get_url = 'http://gene-regulation.com/cgi-bin/pub/programs/patch/bin/patch.cgi'
#! 使用cookie文件进行登录
cookie_filename = os.path.join(os.path.dirname(__file__), 'cookie.txt')
cookie_aff = http.cookiejar.MozillaCookieJar(cookie_filename)
cookie_aff.load(cookie_filename, ignore_discard=True, ignore_expires=True)
handler = urllib.request.HTTPCookieProcessor(cookie_aff)
opener = urllib.request.build_opener(handler)
# 构造请求头，伪装成浏览器
user_agent = r'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_13_3) AppleWebKit/537.36' \
             r' (KHTML, like Gecko) Chrome/61.0.3163.79 Safari/537.36'
headers = {'User-Agent': user_agent, 'Connection': 'keep-alive'}
# 构造post请求post表单
searchvalue = {"Status": "First",
               'searchName': 'default',
               'usr_seq': 'default',
               "seqStat": "DEL",
               'sequenceName': 'default.seq',
               "site_opt": "OUR",
               'group': 'plants',
               'minLen': 8,
               'mismatch': 1,
               'penalty': 100,
               'boundary': 87.5}
# 初始化序列信息
searchvalue['theSequence'] = ''
# 读取基因序列信息，赋值给字典
genelist = fastaread(sys.argv[1])
# 创建输出文件句柄
patchout = open(sys.argv[2], 'w')
# 用于记录序列数目
flag = 0
# 循环添加序列信息
for gene in tqdm(genelist.keys(), desc="request is doing"):
    flag += 1
    # 拼接字符串操作
    searchvalue['theSequence'] = '%s%s%s%s%s%s' % (
        searchvalue['theSequence'], ">", gene, " \n", genelist[gene], "\n")
    # 每隔200个序列发起一次请求，但是最后还会剩下不能够整除200的一些序列
    if(flag % 200 == 0):
        # 对post内容进行编码
        searchtdata = urllib.parse.urlencode(searchvalue).encode()
        # 使用cookie登陆get_url
        get_request = urllib.request.Request(
            get_url, searchtdata, headers=headers)
        get_response = opener.open(get_request)
        # 创建 BeautifulSoup对象
        soup = BeautifulSoup(get_response.read().decode(),
                             features="html.parser")
        # BeautifulSoup 解析说明文档https://www.crummy.com/software/BeautifulSoup/bs4/doc.zh
        for form in soup.find_all("form", action="/cgi-bin/pub/programs/patch/bin/files.cgi"):
            # 查找所有的提交表单
            if(form.previous_sibling.previous_sibling.previous_sibling.string != None):
                patchout.write(
                    ">"+form.previous_sibling.previous_sibling.previous_sibling.string+"\n")
            else:
                print("当前的"+flag+"存在没有的结果！")
            # 根据兄弟节点获得每个序列id和结果信息
            pre = form.next_sibling.next_sibling
            # 如果没有结果，它的string就为no result find
            if(pre.string != None):
                patchout.write(pre.string+"\n")
            else:
                for string in pre.strings:
                    patchout.write(string)
                patchout.write("\n")
        searchvalue['theSequence'] = ''
        # 推迟1s，发送请求
        time.sleep(1)
        # 当只剩下最后89不能够满足整除的条件时
    elif(flag == len(genelist)):
        # 对post内容进行编码
        searchtdata = urllib.parse.urlencode(searchvalue).encode()
        # 使用cookie登陆get_url
        get_request = urllib.request.Request(
            get_url, searchtdata, headers=headers)
        get_response = opener.open(get_request)
        # 创建 BeautifulSoup对象
        soup = BeautifulSoup(get_response.read().decode(),
                             features="html.parser")
        # BeautifulSoup 解析说明文档https://www.crummy.com/software/BeautifulSoup/bs4/doc.zh
        for form in soup.find_all("form", action="/cgi-bin/pub/programs/patch/bin/files.cgi"):
            # 查找所有的提交表单
            if(form.previous_sibling.previous_sibling.previous_sibling.string != None):
                patchout.write(
                    ">"+form.previous_sibling.previous_sibling.previous_sibling.string+"\n")
            else:
                print("当前的"+flag+"没有结果！")
            # 根据兄弟节点获得每个序列id和结果信息
            pre = form.next_sibling.next_sibling
            # 如果没有结果，它的string就为no result find
            if(pre.string != None):
                patchout.write(pre.string+"\n")
            else:
                for string in pre.strings:
                    patchout.write(string)
                patchout.write("\n")
    else:
        continue

patchout.close()

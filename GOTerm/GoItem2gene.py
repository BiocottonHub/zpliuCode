'''/*
 * @Author: zpliu 
 * @Date: 2020-08-07 08:39:47 
 * @Last Modified by:   zpliu 
 * @Last Modified time: 2020-08-07 08:39:47 
 */
'''

import urllib.error
import urllib.request
import urllib.parse
import json
import time
import sys


'''
伪装成浏览器
'''
user_agent = r'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_13_3) AppleWebKit/537.36' \
             r' (KHTML, like Gecko) Chrome/61.0.3163.79 Safari/537.36'
headers = {'User-Agent': user_agent, 'Connection': 'keep-alive'}


def getGoitemBykeyword(keyword,):
    count = '1'
    getGoItemUrl = 'http://golr.geneontology.org/solr/select?defType=edismax&qt=standard&indent=on&wt=json&rows='+count+'&start=0&fl=annotation_class,description,source,idspace,synonym,alternate_id,annotation_class_label,score,id&facet=true&facet.mincount=1&facet.sort=count&json.nl=arrarr&facet.limit=25&hl=true&hl.simple.pre=%3Cem%20class=%22hilite%22%3E&hl.snippets=1000&fq=document_category:%22ontology_class%22&fq=idspace:%22GO%22&fq=is_obsolete:%22false%22&facet.field=source&facet.field=idspace&facet.field=subset&facet.field=is_obsolete&q=' + \
        keyword+'*&qf=annotation_class%5E3&qf=annotation_class_label_searchable%5E5.5&qf=description_searchable%5E1&qf=synonym_searchable%5E1&qf=alternate_id%5E1&packet=1&callback_type=search'
    try:
        response = urllib.request.urlopen(getGoItemUrl)  # 发起get请求
        jsonData = json.loads(response.read().decode('utf-8'))  # 解析成json数据
        # time.sleep(1)  # 缓一缓
        count = str(jsonData['response']['numFound'])
        getGoItemUrl = 'http://golr.geneontology.org/solr/select?defType=edismax&qt=standard&indent=on&wt=json&rows='+count+'&start=0&fl=annotation_class,description,source,idspace,synonym,alternate_id,annotation_class_label,score,id&facet=true&facet.mincount=1&facet.sort=count&json.nl=arrarr&facet.limit=25&hl=true&hl.simple.pre=%3Cem%20class=%22hilite%22%3E&hl.snippets=1000&fq=document_category:%22ontology_class%22&fq=idspace:%22GO%22&fq=is_obsolete:%22false%22&facet.field=source&facet.field=idspace&facet.field=subset&facet.field=is_obsolete&q=' + \
            keyword+'*&qf=annotation_class%5E3&qf=annotation_class_label_searchable%5E5.5&qf=description_searchable%5E1&qf=synonym_searchable%5E1&qf=alternate_id%5E1&packet=1&callback_type=search'
        try:
            response = urllib.request.urlopen(getGoItemUrl)  # 发起get请求
            jsonData = json.loads(response.read().decode('utf-8'))  # 解析成json数据
            # time.sleep(1)  # 缓一缓
            # print(jsonData['response']['docs'])
            return jsonData['response']['docs']  # 得到GO item信息
        except urllib.error.URLError as e:
            print("get请求出错:"+e.reason)
    except urllib.error.URLError as e:
        print("get请求出错:"+e.reason)


def getGeneIdByGoItem(GoItem):
    count = '1'
    getGeneInfoUrl = 'http://golr.geneontology.org/solr/select?defType=edismax&qt=standard&indent=on&wt=json&rows='+count+'&start=0&fl=bioentity,bioentity_name,taxon,panther_family,type,source,synonym,bioentity_label,taxon_label,panther_family_label,score,id&facet=true&facet.mincount=1&facet.sort=count&json.nl=arrarr&facet.limit=25&hl=true&hl.simple.pre=%3Cem%20class=%22hilite%22%3E&hl.snippets=1000&fq=regulates_closure:%22' + \
        GoItem+'%22&fq=document_category:%22bioentity%22&facet.field=source&facet.field=taxon_subset_closure_label&facet.field=type&facet.field=panther_family_label&facet.field=annotation_class_list_label&facet.field=regulates_closure_label&q=*:*&packet=1&callback_type=search'
    try:
        response = urllib.request.urlopen(getGeneInfoUrl)  # 发起get请求
        jsonData = json.loads(response.read().decode('utf-8'))  # 解析成json数据
        # time.sleep(1)  # 缓一缓
        count = str(jsonData['response']['numFound'])
        getGeneInfoUrl = 'http://golr.geneontology.org/solr/select?defType=edismax&qt=standard&indent=on&wt=json&rows='+count+'&start=0&fl=bioentity,bioentity_name,taxon,panther_family,type,source,synonym,bioentity_label,taxon_label,panther_family_label,score,id&facet=true&facet.mincount=1&facet.sort=count&json.nl=arrarr&facet.limit=25&hl=true&hl.simple.pre=%3Cem%20class=%22hilite%22%3E&hl.snippets=1000&fq=regulates_closure:%22' + \
            GoItem+'%22&fq=document_category:%22bioentity%22&facet.field=source&facet.field=taxon_subset_closure_label&facet.field=type&facet.field=panther_family_label&facet.field=annotation_class_list_label&facet.field=regulates_closure_label&q=*:*&packet=1&callback_type=search'
        try:
            response = urllib.request.urlopen(getGeneInfoUrl)  # 发起get请求
            jsonData = json.loads(response.read().decode('utf-8'))  # 解析成json数据
            # time.sleep(1)  # 缓一缓
            # print(jsonData['response']['docs'][0])  # 获取GO info信息
            return jsonData['response']['docs']
        except urllib.error.URLError as e:
            print("get请求出错:"+e.reason)
    except urllib.error.URLError as e:
        print("get请求出错:"+e.reason)


searchKeyWord = sys.argv[1]  # 搜索关键字
GoItem = getGoitemBykeyword(searchKeyWord)  # 获取对应的Go item
out = []
# 获取 GO item的基因注释信息
GoKeys = ['id', 'annotation_class_label', 'description', 'source', 'score']
geneInfoKeys = ['bioentity', 'bioentity_label', 'bioentity_name',
                'source', 'type', 'taxon', 'taxon_label', 'panther_family', 'score']
for item in GoItem:
    print("爬取GO编号:\t"+item['id'])
    info = getGeneIdByGoItem(item['id'])
    GoInfo = "\t".join([str(item[i]) for i in GoKeys])  # 获取Go item 描述信息
    if(len(info) == 0):
        geneMessage = "\t".join(
            ['NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA'])
        out.append(GoInfo+"\t"+geneMessage+"\n")
    else:
        for geneItem in info:
            geneMessage = []
          # 有不存在值的列,就输出NA信息
            for key in geneInfoKeys:
                try:
                    geneMessage.append(geneItem[key])
                except KeyError:
                    geneMessage.append("NA")
            out.append(GoInfo+"\t"+"\t".join([str(i)
                                              for i in geneMessage])+"\n")


with open(sys.argv[2], 'w') as File:
    File.write("\t".join(GoKeys)+"\t"+"\t".join(geneInfoKeys)+"\n")
    for item in out:
        File.write(item)

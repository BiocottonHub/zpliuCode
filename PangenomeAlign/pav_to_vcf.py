'''
Descripttion: 
version: 
Author: zpliu
Date: 2023-02-15 14:18:05
LastEditors: zpliu
LastEditTime: 2023-02-16 09:16:41
@param: 
reference:  https://github.com/wjian8/psvcp_v1.01 
'''
import pandas as pd 
import sys 
chromId=sys.argv[1]

#---------------------------------------------------
#TODO 修改样本名称
#---------------------------------------------------
renameSample=pd.read_csv(
    "/data/cotton/MaojunWang/WMJ_fiberFullPopulationRNAseq/GenomeSNPs/RNA_SampleNames_DNAnames.376.txt",
        header=None,sep="\t",index_col=None
    )
renameSample.columns=['RNAsample','ID','DNAsample']    
renameDict={}
for value in renameSample.values:
    renameDict[value[2]+"."]=renameDict.get(
        value[2]+".",value[0]
    )
    
genotypeMap={
    'C':'0/0',
    'A':'1/1'
}
#* 修改染色体编号
ChromDict={
    'Ghir_A01':'1','Ghir_A02':'2','Ghir_A03':'3','Ghir_A04':'4','Ghir_A05':'5','Ghir_A06':'6',
    'Ghir_A07':'7','Ghir_A08':'8','Ghir_A09':'9','Ghir_A10':'10','Ghir_A11':'11','Ghir_A12':'12',
    'Ghir_A13':'13',
    'Ghir_D01':'14','Ghir_D02':'15','Ghir_D03':'16','Ghir_D04':'17','Ghir_D05':'18','Ghir_D06':'19',
    'Ghir_D07':'20','Ghir_D08':'21','Ghir_D09':'22','Ghir_D10':'23','Ghir_D11':'24','Ghir_D12':'25',
    'Ghir_D13':'26',
}


#------------------------------------------------------
#TODO 修改文件列名，并根据MAF进行过滤
#------------------------------------------------------
pav_mergeData=pd.read_csv(
    "read_count_chr/{}/pav_{}".format(chromId,chromId),header=0,index_col=None,sep="\t"
    )
maf=0.05
#* 样本总数
sampleCount=pav_mergeData.shape[1]-2
#* 根据MAF筛选行
pav_mergeData['maf']=pav_mergeData.apply(lambda x:len([i for i in x if i=='A'])/sampleCount,axis=1)
pav_mergeData=pav_mergeData.loc[pav_mergeData['maf'].ge(maf)&pav_mergeData['maf'].lt(1-maf)]   
#* 删除maf列 
pav_mergeData=pav_mergeData.iloc[:,0:-1]

#* 将AC基因型改为0/1 1/1类型
pav_mergeData=pav_mergeData.applymap(lambda x:genotypeMap.get(x,x)) 


#* 添加信息
pav_mergeData.insert(2,'FORMAT','GT')
pav_mergeData.insert(2,'INFO','INFO')
pav_mergeData.insert(2,'FILTER','PASS')
pav_mergeData.insert(2,'QUAL','NA')
pav_mergeData.insert(2,'ALT','A')
pav_mergeData.insert(2,'REF','C')
pav_mergeData['pos'] = pav_mergeData['pos'].map(lambda x:str(x))
#* 染色体编号替换为数字
pav_mergeData.rename(columns={'chr':'#CHROM'},inplace=True)
pav_mergeData.rename(columns={'pos':'POS'},inplace=True)
pav_mergeData['POS'] = pav_mergeData['POS'].map(lambda x:str(x))
pav_mergeData['#CHROM'] = pav_mergeData['#CHROM'].map(ChromDict)
pav_mergeData.insert(2,'ID',pav_mergeData['#CHROM']+"_"+pav_mergeData['POS']+"_pav")
#* 修改样本对应的编号
pav_mergeData.columns=[renameDict.get(i,i) for i in pav_mergeData.columns] 

#* 保存文件
pav_mergeData.to_csv(
    "./PAV_merge/{}_rnaSample_PAV.vcf".format(chromId),
    header=True,index=False,sep="\t")
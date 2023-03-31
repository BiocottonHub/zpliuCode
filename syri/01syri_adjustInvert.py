'''
Descripttion: 
version: 
Author: zpliu
Date: 2023-03-30 17:21:34
LastEditors: zpliu
LastEditTime: 2023-03-31 19:27:06
@param: 
'''
#! 由于syri在inversion数目比较多的时候，会出现报错，为了能够正常的跑；
#! 将query和reference的链修改
#TODO 因此需要根据invert后的情况，将syri跑的结果中坐标进行调整 
import pandas as pd 
import re 
import pysam
import sys 
#query genome
genomeObject=pysam.FastaFile("/cotton/Liuzhenping/Pan-genome/HC04/HC04_V2/HC04_chr_adjust.fa")
#* syri result
#! syri中的坐标是1-base, `testsyri.out`
syriOut=pd.read_csv(sys.argv[1],header=None,index_col=None,sep="\t",low_memory=False)
syriOut=syriOut.loc[(syriOut[5]!="-")&(syriOut[0]!="-")]
syriOut=syriOut.astype(
    {
        1:int,
        2:int,
        6:int,
        7:int
    }
)

def annotateType(x):
    '''替换原有的SVs类型
    '''
    if re.match('INV',x):
        if re.match('INV[0-9]*$',x):
            return re.sub('INV','SYN',x)
        elif re.match('INVAL',x):
            return re.sub('INV','SYN',x)
        elif re.match('INVTR',x):
            return re.sub('INVTR','TRANS',x)
        elif re.match('INVDP',x):
            return re.sub('INVDP','DUP',x)
    elif re.match('SYN',x):
        return re.sub('SYN','INV',x) 
    elif re.match('TRANS',x):
        return re.sub('TRANS','INVTR',x)  
    elif re.match('DUP',x):
        return re.sub('DUP','INVDP',x)  
    else:
        return x  

#* 根据染色体长度对坐标进行互补
def rever_local(chromosomeLength,site):
    return chromosomeLength-site+1
baseDict={
    'A':'T',
    'T':'A',
    'G':'C',
    'C':'G',
}

def getBaseInfor(baseInfo,SVType,matchArray,start,end,genomeObject):
    #TODO 获取序列信息
    #* inv的话会将序列反过来
        if baseInfo!="-":
            baseInfo=genomeObject.fetch(reference=matchArray[5],start=start-1,end=end)
        #* 判断是否需要互补该基因型
        if  re.match("INV",SVType) and len(baseInfo)==1 and baseInfo!="-":
                baseInfo=baseDict.get(baseInfo)
        return baseInfo



def getLocation(PAF_reversed,ChromosomeReverse,chromosomeLength,AlignmentType,genomeObject,syri_queryStart,syri_queryEnd,paf_queryStart,paf_queryEnd,matchArray):
    if PAF_reversed==True and ChromosomeReverse==True:
        #* query染色体在翻转后，部分染色体的PAF需要修改链才能跑出结果
        #? 首先在染色体翻转的情况下，得到翻转下的坐标；然后再进行反向互补
        if AlignmentType=="parent" :
            site1=rever_local(chromosomeLength,syri_queryStart)
            site2=rever_local(chromosomeLength,syri_queryEnd)
            start=min(site1,site2)
            end=max(site1,site2)
            #! SVs的类型不需要进行修改
            return (
               matchArray[0],matchArray[1],
                matchArray[2],matchArray[3],
                matchArray[4],matchArray[5],
                start,end,
                matchArray[8],matchArray[9],
                matchArray[10],matchArray[11] 
            )
        elif AlignmentType=="Alignment":
            site1=rever_local(chromosomeLength,syri_queryStart)
            site2=rever_local(chromosomeLength,syri_queryEnd)
            start=min(site1,site2)
            end=max(site1,site2) 
            return (
               matchArray[0],matchArray[1],
                matchArray[2],matchArray[3],
                matchArray[4],matchArray[5],
                start,end,
                matchArray[8],matchArray[9],
                matchArray[10],matchArray[11] 
            )
        else:
            #* 首先需要对坐标进行翻转，再翻转后在进行链的互补
            invertStart=-1*(syri_queryStart-(paf_queryStart+1))+paf_queryEnd
            invertEnd=-1*(syri_queryEnd-(paf_queryStart+1))+paf_queryEnd
            site1=rever_local(chromosomeLength,invertStart)
            site2=rever_local(chromosomeLength,invertEnd)
            start=min(site1,site2)
            end=max(site1,site2)
            #* 根据SVs的类型判断第4类的碱基为SNP时是否需要翻转
            SVsType=matchArray[9]
            baseInfo=matchArray[4]
            baseInfo=getBaseInfor(baseInfo,SVsType,matchArray,start,end,genomeObject)
            return (
                matchArray[0],matchArray[1],
                matchArray[2],matchArray[3],
                baseInfo,matchArray[5],
                start,end,
                matchArray[8],matchArray[9],
                matchArray[10],matchArray[11] 
            )
    elif PAF_reversed==False and ChromosomeReverse==True:
        #* query染色体在翻转后，染色体的PAF不需要修改链,
        if AlignmentType=="parent" :
            site1=rever_local(chromosomeLength,syri_queryStart)
            site2=rever_local(chromosomeLength,syri_queryEnd)
            start=min(site1,site2)
            end=max(site1,site2)
            return (
               matchArray[0],matchArray[1],
                matchArray[2],matchArray[3],
                matchArray[4],matchArray[5],
                start,end,
                annotateType(matchArray[8]),annotateType(matchArray[9]),
                annotateType(matchArray[10]),matchArray[11],
            )

        elif AlignmentType=="Alignment":
            site1=rever_local(chromosomeLength,syri_queryStart)
            site2=rever_local(chromosomeLength,syri_queryEnd)
            start=min(site1,site2)
            end=max(site1,site2)
            return (
               matchArray[0],matchArray[1],
                matchArray[2],matchArray[3],
                matchArray[4],matchArray[5],
                start,end,
                annotateType(matchArray[8]),annotateType(matchArray[9]),
                annotateType(matchArray[10]),matchArray[11],
            )
        else:
            #* 首先需要对坐标进行翻转，再翻转后在进行链的互补
            site1=rever_local(chromosomeLength,syri_queryStart)
            site2=rever_local(chromosomeLength,syri_queryEnd)
            start=min(site1,site2)
            end=max(site1,site2)
            #* 矫正后的SVs类型
            SVsType=annotateType(matchArray[9])
            baseInfo=matchArray[4]
            baseInfo=getBaseInfor(baseInfo,SVsType,matchArray,start,end,genomeObject)
            return (
               matchArray[0],matchArray[1],
                matchArray[2],matchArray[3],
                baseInfo,matchArray[5],
                start,end,
                annotateType(matchArray[8]),annotateType(matchArray[9]),
                annotateType(matchArray[10]),matchArray[11],
            )
    elif PAF_reversed==True and ChromosomeReverse==False:
        #* query 染色体不经过翻转，但是其PAF需要修改链才能继续跑出结果
        if AlignmentType=="parent" :
            start=syri_queryStart
            end=syri_queryEnd 
            return (
               matchArray[0],matchArray[1],
                matchArray[2],matchArray[3],
                matchArray[4],matchArray[5],
                start,end,
                annotateType(matchArray[8]),annotateType(matchArray[9]),
                annotateType(matchArray[10]),matchArray[11],
            )
        elif AlignmentType=="Alignment":
            start=syri_queryStart
            end=syri_queryEnd
            return (
               matchArray[0],matchArray[1],
                matchArray[2],matchArray[3],
                matchArray[4],matchArray[5],
                start,end,
                annotateType(matchArray[8]),annotateType(matchArray[9]),
                annotateType(matchArray[10]),matchArray[11],
            )
        else:
            #* 首先需要对坐标进行翻转
            invertStart=-1*(syri_queryStart-(paf_queryStart+1))+paf_queryEnd
            invertEnd=-1*(syri_queryEnd-(paf_queryStart+1))+paf_queryEnd
            start=min(invertStart,invertEnd)
            end=max(invertStart,invertEnd)
            #* 矫正后的SVs类型
            SVsType=annotateType(matchArray[9])
            baseInfo=matchArray[4]
            baseInfo=getBaseInfor(baseInfo,SVsType,matchArray,start,end,genomeObject)
            return (
               matchArray[0],matchArray[1],
                matchArray[2],matchArray[3],
                baseInfo,matchArray[5],
                start,end,
                annotateType(matchArray[8]),annotateType(matchArray[9]),
                annotateType(matchArray[10]),matchArray[11],
            )
    else:
        #* 正常的情况, 不用使用该脚本处理
        pass 

#TODO 处理每个共线性区块


out=[]
#----------------------------------------------
#! 染色体比对时是否翻转了
# chromosomeLength=119507322
chromosomeLength=int(sys.argv[2])
if chromosomeLength>0:
    ChromosomeReverse=True
else:
    ChromosomeReverse=False
#! PAF的链情况是否翻转了
PAF_reversed=sys.argv[3]
if PAF_reversed=='1':
    PAF_reversed=True
else:
    PAF_reversed=False
#----------------------------------------------

#* syri的结果是按照顺序排列的
paf_queryStart=0
paf_queryEnd=0
for matchArray in syriOut.values:
    if matchArray[3]=='-' and matchArray[4]=="-" and matchArray[9]!="-":
        #* 该列为Alignment不需要修改对应的坐标
        paf_queryStart=min(int(matchArray[6]),int(matchArray[7]))-1
        paf_queryEnd=max(int(matchArray[6]),int(matchArray[7])) 
        newSyriStart=min(int(matchArray[6]),int(matchArray[7]))
        newSyriEnd=max(int(matchArray[6]),int(matchArray[7]))
        result=getLocation(
            PAF_reversed=PAF_reversed,
            ChromosomeReverse=ChromosomeReverse,
            chromosomeLength=chromosomeLength,
            AlignmentType='Alignment',
            genomeObject=genomeObject,
            syri_queryStart=newSyriStart,
            syri_queryEnd=newSyriEnd,
            paf_queryStart=paf_queryStart,
            paf_queryEnd=paf_queryEnd,
            matchArray=matchArray
        )
    elif matchArray[3]=='-' and matchArray[4]=="-" and matchArray[9]=="-":
        #*该列为parent信息列
        newSyriStart=min(int(matchArray[6]),int(matchArray[7]))
        newSyriEnd=max(int(matchArray[6]),int(matchArray[7])) 
        result=getLocation(
            PAF_reversed=PAF_reversed,
            ChromosomeReverse=ChromosomeReverse,
            chromosomeLength=chromosomeLength,
            AlignmentType='parent',
            genomeObject=genomeObject,
            syri_queryStart=newSyriStart,
            syri_queryEnd=newSyriEnd,
            paf_queryStart=paf_queryStart,
            paf_queryEnd=paf_queryEnd,
            matchArray=matchArray
        )
    else:
        #* 大区块内的小片段
        syri_queryStart=int(matchArray[6])
        syri_queryEnd=int(matchArray[7])
        result=getLocation(
            PAF_reversed=PAF_reversed,
            ChromosomeReverse=ChromosomeReverse,
            chromosomeLength=chromosomeLength,
            AlignmentType='item',
            genomeObject=genomeObject,
            syri_queryStart=syri_queryStart,
            syri_queryEnd=syri_queryEnd,
            paf_queryStart=paf_queryStart,
            paf_queryEnd=paf_queryEnd,
            matchArray=matchArray
        )
    #todo 修改对应的SVs类型：
    out.append(
            result
        )
out=pd.DataFrame(out)


#* 存在部分重叠的情况，导致SNP有问题；过滤掉
out=out.loc[~((out[3]==out[4])&(out[3]!="-"))]
#* 最终的结果
out.to_csv(sys.argv[4],header=False,index=False,sep="\t")
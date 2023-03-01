'''
Descripttion: 
version: 
Author: zpliu
Date: 2023-02-28 22:24:34
LastEditors: zpliu
LastEditTime: 2023-02-28 22:38:57
@param: 
'''
from tempfile import NamedTemporaryFile
import pandas as pd
import pysam
import pybedtools
import sys


# * 每个转座子的数据
class TE_aligment(object):
    def __init__(self, TEid):
        # * pysam.AlignedSegment
        self.TE_UP = []
        self.TE_DOWN = []

    def set_TE_UP(self, AligmentObject):
        self.TE_UP.append(
            AligmentObject
        )

    def set_TE_DOWN(self, AligmentObject):
        self.TE_DOWN.append(
            AligmentObject
        )

    def fetch_TE_UP(self):
        for i in self.TE_UP:
            yield i

    def fetch_TE_DOWN(self):
        for i in self.TE_DOWN:
            yield i


InputPath = sys.argv[1]
sampleName = sys.argv[2]
outPath = sys.argv[3]
# TODO
# * 区分转座子，并统计每个坐标上mapping的reads数目
DownbamFile = pysam.AlignmentFile(
    "{}/{}_TE-down.bam".format(InputPath, sampleName), 'rb'
)
UpbamFile = pysam.AlignmentFile(
    "{}/{}_TE-up.bam".format(InputPath, sampleName), 'rb'
)

# * 储存TE的比对信息
TE_objectDict = {}
for i in DownbamFile.fetch():
    TEid = i.qname.split(":")[0]
    TE_objectDict[TEid] = TE_objectDict.get(TEid,
                                            TE_aligment(TEid)
                                            )
    TE_objectDict[TEid].set_TE_DOWN(i)

for i in UpbamFile.fetch():
    TEid = i.qname.split(":")[0]
    TE_objectDict[TEid] = TE_objectDict.get(TEid,
                                            TE_aligment(TEid)
                                            )
    TE_objectDict[TEid].set_TE_UP(i)

# * 临时文件
tmp_TE_UP_BAM = NamedTemporaryFile(mode='w+')
tmp_TE_Down_BAM = NamedTemporaryFile(mode='w+')
tmp_TE_UP_Bed = NamedTemporaryFile(mode='w+')
tmp_TE_Down_Bed = NamedTemporaryFile(mode='w+')

def site_readCover(BAMfile,outBedFile):
    '''
    #TODO 关于每个位点被几个read比对上的问题
    #! 统计每个碱基read count
    '''
    #---------------------------------------
    bamFile=pysam.AlignmentFile(BAMfile,'rb',require_index=False)
    siteReadCount=[]
    for i in bamFile.pileup():
            if i.n<60:
                siteReadCount.append(
                    ( i.reference_name,i.pos+1,i.pos+1,i.n)
                )
    siteReadCount=pd.DataFrame(siteReadCount)
    if siteReadCount.shape[0]==0:
        #* TE reads mapping情况不满足count < 60
        return False 
    siteReadCount.sort_values(by=[0,1])
    #* read信息统计文件
    siteReadCount.to_csv(outBedFile,header=False,index=False,sep="\t") 
    return True

#* 所有TE的数据
TE_insertData=[]
for TEid in TE_objectDict.keys():
    print("deal with TE {}".format(TEid))
    with pysam.AlignmentFile(tmp_TE_Down_BAM.name, "wb", header=DownbamFile.header) as outf:
        #* 将该转座子的所有信息保存到BAM文件中
        for i in TE_objectDict[TEid].fetch_TE_DOWN():
            outf.write(i)
    with pysam.AlignmentFile(tmp_TE_UP_BAM.name, "wb", header=UpbamFile.header) as outf:
        #* 将该转座子的所有信息保存到BAM文件中
        for i in TE_objectDict[TEid].fetch_TE_UP():
            outf.write(i)
    #* 获取BAM文件中每个坐标的reads覆盖度
    if not site_readCover(tmp_TE_Down_BAM.name,tmp_TE_Down_Bed.name) or not site_readCover(tmp_TE_UP_BAM.name,tmp_TE_UP_Bed.name):
        ## reads 数目太少，无法继续进行merge
        tmp_TE_Down_BAM.truncate()
        tmp_TE_UP_BAM.truncate()
        tmp_TE_Down_Bed.truncate()
        tmp_TE_UP_Bed.truncate() 
        continue
    DownBed=pybedtools.BedTool(tmp_TE_Down_Bed.name)
    UPBed=pybedtools.BedTool(tmp_TE_UP_Bed.name)
    #* 在连续的位置上有插入
    DownMerge=DownBed.merge(c=4,o='max')
    UPMerge=UPBed.merge(c=4,o='max')
    InterSectBed=UPMerge.intersect(
        DownMerge,wo=True
    )
    for i in InterSectBed:
        TE_insertData.append(
            [TEid]+str(i).strip("\n").split("\t")
        )
    tmp_TE_Down_BAM.truncate()
    tmp_TE_UP_BAM.truncate()
    tmp_TE_Down_Bed.truncate()
    tmp_TE_UP_Bed.truncate() 

TE_insertData=pd.DataFrame(TE_insertData)
TE_insertData=TE_insertData.astype(
    {
        2:int,
        3:int,
        4:int,
        6:int,
        7:int,
        8:int,
        9:int
    }
)
TE_insertData.to_csv(
    "{}/{}-intersectReadCount.txt".format(
        outPath,sampleName
    ),
    header=False,index=False,sep="\t"
)
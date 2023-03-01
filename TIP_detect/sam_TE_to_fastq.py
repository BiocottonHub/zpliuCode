'''
Descripttion: 
version: 
Author: zpliu
Date: 2023-02-27 20:57:52
LastEditors: zpliu
LastEditTime: 2023-03-01 18:19:06
@param: 
'''
import sys 
import pysam 
#* 最终保存文件
#? J100_TE.fastq
outFile='{}/{}_TE-split.fastq'.format(sys.argv[2],sys.argv[3])
samfile = pysam.AlignmentFile("-", "r")
with open(outFile,'w') as FASTQ:
    for i in samfile:
        readName=i.qname
        refName=i.reference_name
        readSeq=i.query_sequence
        seqQuality=i.to_string().split("\t")[10]
        #* 修改read name
        readName=readName.split(":")
        readName[0]=refName
        readName=":".join(readName)
        FASTQ.write(
            "@{}\n{}\n+\n{}\n".format(
                readName,readSeq,seqQuality
            )
        )
###
 # @Descripttion: 
 # @version: 
 # @Author: zpliu
 # @Date: 2023-03-01 18:05:57
 # @LastEditors: zpliu
 # @LastEditTime: 2023-03-01 18:11:49
 # @@param: 
### 
#TODO 根据split-reads merge后的BAM文件进行统计
#! 每个样本大概耗时25分钟
inputPath='/public/home/zpliu/TIP/Detect_TIPs/split_reads_mapping/test/'
outPath='/public/home/zpliu/TIP/Detect_TIPs/split_reads_mapping/test/'
sampleName='ZY45'
python readCluster_to_TE.py  ${inputPath} ${sampleName} ${outPath} 

#* 删除文件？

rm  ${inputPath}/${sampleName}_TE-down.bam* -rf 
rm  ${inputPath}/${sampleName}_TE-up.bam* -rf 

###
 # @Descripttion: 
 # @version: 
 # @Author: zpliu
 # @Date: 2023-02-28 20:46:54
 # @LastEditors: zpliu
 # @LastEditTime: 2023-03-01 18:20:47
 # @@param: 
### 
sampleName='ZY45'
inputPath='./'
outPath='./'

#*TODO 合并所有的BAM文件
samtools merge -f -u ${outPath}/${sampleName}_TE-split-5.bam ${inputPath}/${sampleName}_TE-splitjunction-5-*
samtools merge -f -u ${outPath}/${sampleName}_TE-split-3.bam ${inputPath}/${sampleName}_TE-splitjunction-3-* 
#* 按照比对位置进行排序 
samtools sort ${outPath}/${sampleName}_TE-split-3.bam -o ${outPath}/${sampleName}_TE-split-3.bam
samtools sort ${outPath}/${sampleName}_TE-split-5.bam -o ${outPath}/${sampleName}_TE-split-5.bam
# rm  ${sampleName}_TE-splitjunction-[35]-*.bam -rf 


#TODO 提取reads在断点两端的比对情况
samtools view -F 16 -bh ${outPath}/${sampleName}_TE-split-5.bam >${outPath}/${sampleName}_TE-split-5+.bam
samtools view -f 16 -bh ${outPath}/${sampleName}_TE-split-5.bam >${outPath}/${sampleName}_TE-split-5-.bam

samtools view -F 16 -bh ${outPath}/${sampleName}_TE-split-3.bam >${outPath}/${sampleName}_TE-split-3+.bam
samtools view -f 16 -bh ${outPath}/${sampleName}_TE-split-3.bam >${outPath}/${sampleName}_TE-split-3-.bam

#* TE在转座子的结束位置
samtools merge -f -u ${outPath}/${sampleName}_TE-up.bam  ${outPath}/${sampleName}_TE-split-5-.bam ${outPath}/${sampleName}_TE-split-3+.bam
samtools merge -f -u ${outPath}/${sampleName}_TE-down.bam  ${outPath}/${sampleName}_TE-split-5+.bam ${outPath}/${sampleName}_TE-split-3-.bam 
samtools sort ${outPath}/${sampleName}_TE-up.bam -o ${outPath}/${sampleName}_TE-up.bam
samtools sort ${outPath}/${sampleName}_TE-down.bam -o ${outPath}/${sampleName}_TE-down.bam 


samtools index ${outPath}/${sampleName}_TE-down.bam
samtools index ${outPath}/${sampleName}_TE-up.bam 

#* 删除多余文件
rm ${outPath}/${sampleName}_TE-split-5-.bam ${outPath}/${sampleName}_TE-split-5+.bam -rf 
rm ${outPath}/${sampleName}_TE-split-3-.bam ${outPath}/${sampleName}_TE-split-3+.bam -rf  

rm ${inputPath}/${sampleName}_TE-splitjunction-5-* -rf 
rm ${inputPath}/${sampleName}_TE-splitjunction-3-* -rf 



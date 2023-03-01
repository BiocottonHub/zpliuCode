###
 # @Descripttion: 
 # @version: 
 # @Author: zpliu
 # @Date: 2023-03-01 18:21:40
 # @LastEditors: zpliu
 # @LastEditTime: 2023-03-01 18:23:55
 # @@param: 
### 
module load Bowtie2/2.4.1
module load picard/2.23.9

#* 重测序数据路径
rawBase=/cotton/XinyuanChen/NewPhytol/251AccessionsDNAseqData/
#* 重测序编号文件
sampleList=/cotton/Liuzhenping/Pan-genome/Genotype/sampleID_250.txt
#! 将重测序数据Mapping到参考基因组, 并去除PCR重复
bowtieIndex='/public/home/zpliu/TIP/TE_annotion/Ghirsutum.chr.bowtie'
sampleName='ZY446'
fastq1=$rawBase/${sampleName}_1.fq.gz
fastq2=$rawBase/${sampleName}_2.fq.gz
outPath=''

threads=4
bowtie2 --mp 13 --rdg 8,5 --rfg 8,5 --very-sensitive -x ${bowtieIndex} -1 ${fastq1} -2 ${fastq2} -p ${threads} \
        |samtools view -bhS -o ${outPath}/${sampleName}.bam
#* 对输出后BAM进行排序
java -jar $EBROOTPICARD/picard.jar SortSam I=${outPath}/${sampleName}.bam \
            O=${outPath}/${sampleName}_sorted.bam SORT_ORDER=coordinate
#! 对输出后的文件进行PCR去重
java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
        I=${outPath}/${sampleName}_sorted.bam O=${outPath}/${sampleName}_sorted_dedup.bam \
        M=${outPath}/${sampleName}_dedup.metrics REMOVE_DUPLICATES=true

rm  ${outPath}/${sampleName}_sorted.bam ${outPath}/${sampleName}_dedup.metrics ${outPath}/${sampleName}.bam -rf 
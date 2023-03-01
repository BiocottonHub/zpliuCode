###
# @Descripttion:
# @version:
# @Author: zpliu
# @Date: 2023-03-01 18:12:20
 # @LastEditors: zpliu
 # @LastEditTime: 2023-03-01 18:20:31
# @@param:
###
#TODO 提取未比对到参考基因组上的序列
module load Bowtie2/2.4.1
samplName=''
inputPath=''
outPath=''
TE_index='/public/home/zpliu/TIP/Detect_TIPs/TE_indexs/TE_EDTA_intact'
#* 是否多进程
threadsNum=2

samtools view -f 4 -u ${inputPath}/${samplName}_sorted_dedup.bam >${outPath}/${samplName}_filter.bam
java -jar ${EBROOTPICARD}/picard.jar \
    SamToFastq I=${outPath}/${samplName}_filter.bam \
    FASTQ=${outPath}/${samplName}_filter.1.fastq \
    SECOND_END_FASTQ=${outPath}/${samplName}_filter.2.fastq

#* 将双端数据还需要将1和2合并为一个文件
cat ${outPath}/${samplName}_filter.1.fastq ${outPath}/${samplName}_filter.2.fastq >${outPath}/${samplName}_unmapping.fastq
rm ${outPath}/${samplName}_filter.1.fastq ${outPath}/${samplName}_filter.2.fastq -f
rm ${outPath}/${samplName}_filter.bam -rf

bowtie2 -x ${TE_index} \
    -U ${outPath}/${samplName}_unmapping.fastq \
    --local --very-sensitive \
    --threads ${threadsNum} \
    |samtools view -h -F 4 - \
    |awk '$1~/^@/{print $0}$6~/^[2-8][0-9]S/ || $6~/[2-8][0-9]S$/ {print $0}' \
    |python sam_TE_to_fastq.py - ${outPath} ${samplName} 

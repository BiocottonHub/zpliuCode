###
# @Descripttion:
# @version:
# @Author: zpliu
# @Date: 2023-02-28 09:37:46
 # @LastEditors: zpliu
 # @LastEditTime: 2023-03-01 15:52:50
# @@param:
###
#*5.进行迭代mapping
LENGTH=100
length=$(($((LENGTH / 2)) - 1))
referenceIndex='/public/home/zpliu/TIP/TE_annotion/Ghirsutum.chr.bowtie'
inputPath='/public/home/zpliu/TIP/Detect_TIPs/split_reads_mapping/'
sampleName='ZY45'
outPath='./'

module load Bowtie2/2.4.1
#* 第一次5' mapping
#! 这一步耗时5分钟
bowtie2 -x ${referenceIndex} -U ${inputPath}/${sampleName}_TE-split.fastq \
    --un ${outPath}/${sampleName}_TE-split-5-$length -5 $length \
    --mp 13 --rdg 8,5 --rfg 8,5 --very-sensitive --quiet |
    samtools view -Sbu -F 4 |
    samtools sort - -o ${outPath}/${sampleName}_TE-splitjunction-5-$length.bam

#* 这一步大概1个小时？
for ((i = $((length + 1)); $i < $((LENGTH - 18)); i = $i + 1)); do
    echo -n "|"
    previous=$(($i - 1))
    bowtie2 -x ${referenceIndex} -U ${outPath}/${sampleName}_TE-split-5-$previous \
        --un ${outPath}/${sampleName}_TE-split-5-$i \
        -5 $i --mp 13 --rdg 8,5 --rfg 8,5 \
        --very-sensitive --quiet |
        samtools view -Sbu -F 4 |
        samtools sort - -o ${outPath}/${sampleName}_TE-splitjunction-5-$i.bam

    #* 删除之前未mapping上的reads
    rm -f ${sampleName}_TE-split-5-$previous
done
#* 最后一轮驯化剩下的unmapping的修剪reads
rm ${outPath}/${sampleName}_TE-split-5-$(($i - 1)) -rf
echo -n "]"
echo -e "\n"

#* 3'端mapping
length=$(($((LENGTH / 2)) - 1))
bowtie2 -x ${referenceIndex} -U ${inputPath}/${sampleName}_TE-split.fastq \
    --un ${outPath}/${sampleName}_TE-split-3-$length -3 $length \
    --mp 13 --rdg 8,5 --rfg 8,5 --very-sensitive --quiet |
    samtools view -Sbu -F 4 |
    samtools sort - -o ${outPath}/${sampleName}_TE-splitjunction-3-$length.bam

for ((i = $((length + 1)); $i < $((LENGTH - 18)); i = $i + 1)); do
    echo -n "|"
    previous=$(($i - 1))
    bowtie2 -x ${referenceIndex} -U ${outPath}/${sampleName}_TE-split-3-$previous \
        --un ${outPath}/${sampleName}_TE-split-3-$i \
        -3 $i --mp 13 --rdg 8,5 --rfg 8,5 \
        --very-sensitive --quiet |
        samtools view -Sbu -F 4 |
        samtools sort - -o ${outPath}/${sampleName}_TE-splitjunction-3-$i.bam

    #* 删除之前未mapping上的reads
    rm -f ${outPath}/${sampleName}_TE-split-3-$previous
done
rm ${outPath}/${sampleName}_TE-split-3-$(($i - 1)) -rf
echo -n "]"
echo -e "\n"

#todo 删除所有unmapping到参考基因组的soft-colipped read文件
# rm ${sampleName}_TE1.split.fastq -rf



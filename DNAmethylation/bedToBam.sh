###
# @Descripttion:
# @version:
# @Author: zpliu
# @Date: 2021-01-14 10:38:27
# @LastEditors: zpliu
# @LastEditTime: 2021-01-14 11:41:26
# @@param:
###

chromosomeBed='~/github/zpliuCode/circos/A2/chromosome.bed'
bedToolsPath=' ~/software/bedtools2-2.29.0/bin/'
for sample in CpG CHH CHG; do
  bsub -q high -n 1 -R span[hosts=1] -e bedTobam.err -J bedToBam "
  awk '\$1~/Chr/{OFS=\"\t\";printf \$1,\$2,\$3,\$1\"-\"\$2\"-\"\$4\"-\"\$5\"-\"; printf(\"%f\",\$7)}' ${sample}_fdr.bed |${bedToolsPath}/bedToBam -g ${chromosomeBed} -i - | samtools sort -O BAM - >${sample}_fdr_sorted.bam
  samtools index ${sample}_fdr_sorted.bam"

done

for i in 1; do
  bsub -q high -n 1 -R span[hosts=1] -e bedTobam.err -J bedToBam " 
  allsiteFile=$(ls ./ | grep deduplicated.bismark)
  # gunzip ${allsiteFile}
  # if [ $? -eq 0 ]; then
  #   allsiteFile=$(echo ${allsiteFile} | sed 's/\.gz$//')
  # fi
  awk '\$1~/Chr/{OFS=\"\t\";print \$1,\$2,\$3,\$1\"-\"\$2\"-\"\$5\"-\"\$6\"-\"\$4}' ${allsiteFile} |${bedToolsPath}/bedToBam -g ${chromosomeBed} -i - | samtools sort -O BAM - >all_site.bam
  samtools index all_site.bam
"
done

###
# @Descripttion:
# @version:
# @Author: zpliu
# @Date: 2020-12-08 20:17:49
# @LastEditors: zpliu
# @LastEditTime: 2020-12-08 20:47:54
# @@param:
###
Replicate=1
for count in 2; do
  for sample in CHH CHG CpG; do
    bsub -q smp -J splitFile -e bismark.err -o bismark.out -n 1 -R span[hosts=1] "python /data/cotton/zhenpingliu/QingxinSong_GB_DNAmethlation/stat_methlation/splitFile.py tmp/${sample}  ${sample}${count}_ |xargs -I {} bash ./bionormal.sh {}"
  done
done

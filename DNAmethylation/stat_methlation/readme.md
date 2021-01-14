

### 对甲基化read进行统计测验



### 1.统计每个位点的read数目

需要将`bionormal_LSF.sh`和`bionormal.sh`脚本复制到甲基化数据所在目录

+ 其中`bionormal_LSF.sh`是用来对原始文件进行拆分，而`bionormal.sh`根据拆分后的文件，分别对每个文件中位点的read数目进行统计，在新生产的目录下对每个拆分的文件产生`_out`后缀的文件
+ 可以将拆分后的`_out`文件进行合并,可以根据合并后文件的大小判断是否进行第二次拆分
+ 将最后一次拆分得到的文件合并后，进行二项分布测验

```python
##对原始的read文件进行拆分 bionormal_LSF.sh
Replicate=1
for count in 1; do
  for sample in CHH; do
    bsub -q smp -J splitFile -e bismark.err -o bismark.out -n 1 -R span[hosts=1] "python /data/cotton/zhenpingliu/QingxinSong_GB_DNAmethlation/stat_methlation/splitFile.py  ./CHH_context_Rep1.deduplicated.txt.gz  ${sample}${count}_ |xargs -I {} bash ./bionormal.sh {}"
  done
done
##对拆分后的文件进行read统计 bionormal.sh
for i in `ls ${1}`
do
bsub -q high -n 1 -J ${i} -R span[hosts=1] " python /data/cotton/zhenpingliu/QingxinSong_GB_DNAmethlation/stat_methlation/extract_read_count.py ${1}/${i} ${1}/${i}_out"
done

```

#### 1.1将每个单独统计read的文件，进行合并后在统计一次read数目

```bash
##合并所以read count的文件
cat 单独文件* >总文件 
```

### 2.进行二项分布的时候，既可以选择用总的文件进行二项分布也可以将文件拆分后在进行二项分布测验

>  这里将文件进行拆分后，单独对每个文件进行二项分布检验

```bash
##拆分文件的任务 bionormal_LSF.sh
Replicate=1
for count in 2; do
  for sample in CHH CHG CpG; do
    bsub -q smp -J splitFile -e bismark.err -o bismark.out -n 1 -R span[hosts=1] "python /data/cotton/zhenpingliu/QingxinSong_GB_DNAmethlation/stat_methlation/splitFile.py tmp/${sample}  ${sample}${count}_ |xargs -I {} bash ./bionormal.sh {}"
  done
done

##进行二项分布的任务 bionormal.sh

for i in `ls ${1}`
do
bsub -q high -n 1 -J ${i} -R span[hosts=1] " python /data/cotton/zhenpingliu/QingxinSong_GB_DNAmethlation/stat_methlation/bionormal_site.py  -process 1 -input  ${1}/${i} -output  ${1}/${i}_out"
done

```






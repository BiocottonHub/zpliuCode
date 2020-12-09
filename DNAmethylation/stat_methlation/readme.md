

### 对甲基化read进行统计测验

需要将`bionormal_LSF.sh`和`bionormal.sh`脚本复制到甲基化数据所在目录

+ 其中`bionormal_LSF.sh`是用来对原始文件进行拆分，而`bionormal.sh`根据拆分后的文件，分别对每个文件中位点的read数目进行统计，在新生产的目录下对每个拆分的文件产生`_out`后缀的文件
+ 可以将拆分后的`_out`文件进行合并,可以根据合并后文件的大小判断是否进行第二次拆分
+ 将最后一次拆分得到的文件合并后，进行二项分布测验

```python
##
python bionomal_site.py -input 最后一次拆分文件 -process 进程数 -output 输出文件 
```

进行文件拆分后，可以随机提交每个拆分后文件的任务

```bash
Replicate=1
for count in 2; do
  for sample in CHH CHG CpG; do
    bsub -q smp -J splitFile -e bismark.err -o bismark.out -n 1 -R span[hosts=1] "python /data/cotton/zhenpingliu/QingxinSong_GB_DNAmethlation/stat_methlation/splitFile.py tmp/${sample}  ${sample}${count}_ |xargs -I {} bash ./bionormal.sh {}"
  done
done
```






<!--
 * @Descripttion: 
 * @version: 
 * @Author: zpliu
 * @Date: 2023-03-31 11:39:40
 * @LastEditors: zpliu
 * @LastEditTime: 2023-03-31 20:08:20
 * @@param: 
-->
## Syri鉴定SVs时存在的四种情况
1. 染色体没有翻转，同源染色体的syri也能跑出结果
2. 染色体翻转后，同源然染色体syri才能跑出结果
3. 染色体翻转后，同源染色体能跑出结果，但是非同源染色体的比对无法跑出结果


### 1.minimap将基因组序列进行比对
1条染色体和另外一个基因组的所有染色体进行比对

```bash
    minimap2 -t 10 -ax asm5 --eqx ${reference} \
        /public/home/zpliu/Pan-genome/parallele/syri_SVs/HC04_At_chr_split/${query} >${chrom}/${chrom}.sam
    "
```

### 2.将比对的SAM文件转为PAF格式
```bash
export PATH=/cotton/Liuzhenping/Pan-genome/software/minimap2-2.23_x64-linux:$PATH
for query in $(ls ../HC04_At_chr_split); do
    chrom=$(echo ${query} | sed 's/.fa//g')
    paftools.js sam2paf ${chrom}/${chrom}.sam >${chrom}/${chrom}.paf
done
```


### 3.区分reference 染色体进行SVs的鉴定
其中`HC04_${queryChrom}-Allchr_revQuery.paf` 为`1-vs-many` 的paf结果

```bash
chrom=$1
queryChrom=$2
awk -v chr=Chr$chrom '$6==chr{print $0}' ../HC04_${queryChrom}-Allchr_revQuery.paf >HC04_${queryChrom}-Chr${chrom}_revQuery.paf
syri -c HC04_${queryChrom}-Chr${chrom}_revQuery.paf \
    -r /public/home/zpliu/Pan-genome/Genotype_diploid/bwa_index/J85.fa \
    -q /public/home/zpliu/Pan-genome/parallele/syri_SVs/HC04_At_chr_split/HC04_${queryChrom}_rev.fa \
    -F P \
    --prefix HC04_${queryChrom}-Chr${chrom}_revQuery

```

### 4.针对单个reference chromosome进行SVs的鉴定

#### 4.1 当syri处理paf能够正常跑时

`*syri.out`即为最终的结果

```bash
syri -c HC04_${queryChrom}-Chr${chrom}.paf \
            -r /public/home/zpliu/Pan-genome/Genotype_diploid/bwa_index/J85.fa \
            -q /public/home/zpliu/Pan-genome/parallele/syri_SVs/HC04_At_chr_split/HC04_${queryChrom}.fa \
            -F P \
            --prefix syri_HC04_${queryChrom}-Chr${chrom}
```


#### 4.2 当同源染色体需要翻转后才能正常跑syri
例如当前文件`HC04_A01-Chr01_revQuery.paf` 就是A01翻转后与Chr01比对得到的paf文件才能跑syri; 需要将syri处理`HC04_A01-Chr01_revQuery.paf`得到的`HC04_A01-Chr01_revQuerysyri.out` 坐标进行调整

+ 第一个参数是需要修改的`*syri.out` 文件
+ 第二个参数是染色体长度，如果query染色体经过翻转的就设置为真实的长度；否则设置为0
+ 第三个参数是PAF得到`*syri.out`时，对应的PAF是否经过修改；如果修改了就是1，否则为0
+ 第四个参数则是最终输出的`*syri.out` 文件, 该文件是以同源染色体未经过翻转状态下，定义的SVs结果

在这个例子中`HC04_A01-Chr01_revQuery.paf` 文件是可以直接跑syri得到对应的`*syri.out` 文件，所有PAF是未经过修改的

```bash
#* 跑syri的query genome同样使用翻转后的染色体
#* 在矫正syri输出结果时，用参考基因组中的染色体
python  01syri_adjustInvert.py  HC04_A01-Chr01_revQuerysyri.out 119507322 0 HC04_A01-Chr01syri.out
```


#### 4.3 当同源染色体需要翻转后，其他染色体需要修改PAF文件才能正常跑syri

```bash
# todo修改PAF文件
less HC04_A01-Chr03_revQuery.paf |
    awk '$5=="-"{OFS="\t";
    print $1,$2,$3,$4,"+",$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17
        }$5=="+"{OFS="\t";
        print $1,$2,$3,$4,"-",$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17}' >HC04_A01-Chr03_revQuery_changeStand.paf
#* 跑**syri**
syri -c HC04_A01-Chr03_revQuery_changeStand.paf \
    -r /public/home/zpliu/Pan-genome/Genotype_diploid/bwa_index/J85.fa \
    -q /public/home/zpliu/Pan-genome/parallele/syri_SVs/HC04_At_chr_split/HC04_${queryChrom}_rev.fa \
    -F P \
    --prefix HC04_A01-Chr03_revQuery_changeStand

#TODO 根据修改PAF后得到的syri.out 文件进行修改
python  01syri_adjustInvert.py  HC04_A01-Chr03_revQuery_changeStandsyri.out 119507322 1 HC04_A01-Chr03syri.out
```


#### 4.4 当同源染色体不需要翻转，其他染色体需要修改PAF文件才能正常跑syri

```bash
#* 修改PAF文件

cat HC04_A02-Chr03.paf |
    awk '$5=="-"{OFS="\t";
    print $1,$2,$3,$4,"+",$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17
        }$5=="+"{OFS="\t";
        print $1,$2,$3,$4,"-",$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17}' > HC04_A02-Chr03_changeStand.paf
#* 改链后的PAF进行syri鉴定
syri -c HC04_A02-Chr03_changeStand.paf \
    -r /public/home/zpliu/Pan-genome/Genotype_diploid/bwa_index/J85.fa \
    -q /public/home/zpliu/Pan-genome/parallele/syri_SVs/HC04_At_chr_split/HC04_${queryChrom}.fa \
    -F P \
    --prefix HC04_A02-Chr03_changeStand_
python  01syri_adjustInvert.py  HC04_A02-Chr03_changeStand_syri.out 0 1 HC04_A02-Chr03syri.out
```

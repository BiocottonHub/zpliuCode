<!--
 * @Descripttion: 
 * @version: 
 * @Author: zpliu
 * @Date: 2023-03-01 18:29:09
 * @LastEditors: zpliu
 * @LastEditTime: 2023-03-01 18:30:40
 * @@param: 
-->

### 1.aligment sequence to reference

`01-resequence-mapping.sh` 

### 2.extract unmapping reads 

+ extract unmapping reads
+ mapping to TE lib
+ convert bam to fastq

`02-unmapping_reference.sh` 

### 3. iteration mapping potential split-reads to reference 

`03-split.reads-mapping.sh` 

### 4. merge split reads

`04-split.read-merge.sh` 

### identify TE insert site

`05-split-read-cluster.sh` 
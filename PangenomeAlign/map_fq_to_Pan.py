'''
Descripttion: 
version: 
Author: zpliu
Date: 2023-01-28 17:00:51
LastEditors: zpliu
LastEditTime: 2023-01-28 17:28:01
@param: 
'''
#!/usr/local/bin/python3  
# -*- coding: utf-8 -*-  
# Usage: python3  $0 -t 9 -fqd fq_dir -r ReferenceFile -br bam_dir
# Author: Jian Wang
# Email: wjian@gdaas.cn

import os, sys, glob, time, argparse, threading, re

my_parser = argparse.ArgumentParser()
my_parser.add_argument('-t', dest='thread',type=int,default='4') 
my_parser.add_argument('-fqd', dest='input_dir',required=True,help='fq.gz files are in the directory') #
my_parser.add_argument('-r', dest='reference',required=True) #
my_parser.add_argument('-fpf', dest='fq_postfix',default='_1.fq.gz',help="Only fq1 file postfix, like _1.fq.gz or .R1.fastq.gz") #
my_parser.add_argument('-br', dest='bam_dir',default='bam_dir',required=True)
my_parser.add_argument('-bpf', dest='bam_postfix',default='.mapQ20.bam')
my_parser.add_argument('-s', dest='sampleprefix',required=True) 
args = my_parser.parse_args()

ref=args.reference
if not os.path.exists(re.sub('fa','fa.sa',ref)):
    print("bwa index the reference")
    os.system('bwa index %s'%(ref))
if not os.path.exists(re.sub('fa','fa.fai',ref)):
    print("faidx the reference")
    os.system('samtools faidx %s'%(ref))
if not os.path.exists(re.sub('fa','dict',ref)):
    os.system('java -jar ${EBROOTPICARD}/picard.jar CreateSequenceDictionary -R %s -O %s'%(ref,re.sub('fa','dict',ref)))

output_bam_dir =args.bam_dir
if os.path.exists(output_bam_dir):
    sys.exit(output_bam_dir+' exist!!! please check the outputdir')
else:
    os.makedirs(output_bam_dir)
TimeTable = open(args.bam_dir+'TimeTable.txt','w')


fq1_postfix = args.fq_postfix
fq2_postfix = args.fq_postfix.replace('1','2')
sample_name1="{}/{}{}".format(args.input_dir,args.sampleprefix,fq1_postfix)
sample_name2="{}/{}{}".format(args.input_dir,args.sampleprefix,fq2_postfix)






bam_list = []
all_sample_fq = {} 

def bwa(sample_prefix,fq1,fq2,thread): 
    if not os.path.exists(sample_name1) or not os.path.exists(sample_name2) :
            sys.exit(sample_prefix+"' resequence file don't exist!!!")
    else:
        os.system('bwa mem -M -t {5} -R  "@RG\\tID:{0}\\tSM:{0}\\tPL:illumina\\tLB:lib1\\tPU:unit1" {1} {2} {3} | samtools view -bhS -o {4}/{0}.bam'.format(sample_prefix, ref, fq1,fq2, output_bam_dir,thread))
        time.sleep(0.2)
        print(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time()))+'bwa_programming finished***In threading***'+sample_prefix, file=TimeTable)
        TimeTable.flush()
        os.system('java -jar $EBROOTPICARD/picard.jar SortSam INPUT={0}/{1}.bam OUTPUT={0}/{1}_sorted.bam SORT_ORDER=coordinate TMP_DIR=./tmp'.format(output_bam_dir,sample_prefix))
        time.sleep(0.2)
        print(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time()))+'picard_SortSam programming finished***In threading***'+sample_prefix, file=TimeTable)
        TimeTable.flush()
        os.system('rm {0}/{1}.bam'.format(output_bam_dir,sample_prefix))
        os.system('java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups I={0}/{1}_sorted.bam O={0}/{1}_sorted_add.bam RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM={1} TMP_DIR=./tmp'.format(output_bam_dir,sample_prefix))
        time.sleep(0.2)
        print(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time()))+'picard_AddOrReplaceReadGroups programming finished***In threading***'+sample_prefix, file=TimeTable)
        TimeTable.flush()
        os.system('rm {0}/{1}_sorted.bam'.format(output_bam_dir,sample_prefix))
        os.system('java -jar $EBROOTPICARD/picard.jar MarkDuplicates I={0}/{1}_sorted_add.bam O={0}/{1}_sorted_add_dedup.bam M={0}/{1}_sorted_add_dedup.metrics MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=4000 TMP_DIR=./tmp'.format(output_bam_dir,sample_prefix))
        time.sleep(0.2)
        print(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time()))+'picard_MarkDuplicates programming finished***In threading***'+sample_prefix, file=TimeTable)
        TimeTable.flush()
        os.system('rm {0}/{1}_sorted_add.bam'.format(output_bam_dir,sample_prefix))
        os.system('samtools view -h -q 20 -F 4 -F 256 -Sb  {0}/{1}_sorted_add_dedup.bam >{0}/{1}{2}'.format(output_bam_dir,sample_prefix,args.bam_postfix))
        os.system('samtools index {0}/{1}{2} '.format(output_bam_dir,sample_prefix,args.bam_postfix))
        print(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time()))+'samtools_BuildBamIndex programming finished***In threading***'+sample_prefix, file=TimeTable)
        TimeTable.flush()
        os.system('rm {0}/{1}_sorted_add_dedup.bam'.format(output_bam_dir,sample_prefix))
        os.system('rm {0}/{1}_sorted_add_dedup.metrics'.format(output_bam_dir,sample_prefix))

bwa(args.sampleprefix,sample_name1,sample_name2,args.thread)
TimeTable.close()

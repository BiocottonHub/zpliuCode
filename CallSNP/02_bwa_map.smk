rawBase = '/data/cotton/YizanMa/node210_backup/510_resequencing/Trimmed_data'
bwaIndex = '/public/home/weizhang/Call_SNPS/bwa_index/HC04.chr.fa'
platform='ILLUMINA'
threads=5

#'S276',

rule all:
   	input:
	      expand("/public/home/weizhang/Call_SNPS/bam_V2/{sampleName}/{sampleName}_sorted_q20.bam",sampleName=['S131','S276','S420'])


rule bwa_map:
	input:
		fq1="/data/cotton/YizanMa/node210_backup/510_resequencing/Trimmed_data/{sampleName}_paired_R1.fq.gz",
		fq2="/data/cotton/YizanMa/node210_backup/510_resequencing/Trimmed_data/{sampleName}_paired_R2.fq.gz"
        
	output:
		outBam="/public/home/weizhang/Call_SNPS/bam_V2/{sampleName}/{sampleName}_sorted_q20.bam"
	threads:5
	resources:
		mem_mb=15000
	shell:
		"""
		module load SAMtools/1.9
		module load BWA/0.7.17
		bwa mem -Y -K 100000000 -R '@RG\\tID:{wildcards.sampleName}\\tSM:{wildcards.sampleName}\\tPL:{platform}' \
	                -t {threads} {bwaIndex} {input.fq1} {input.fq2} \
	                |samtools  view -bS  -q 20 -@ {threads} \
	                |samtools sort -@ {threads}  -o {output.outBam}
		"""





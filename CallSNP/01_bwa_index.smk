#bwaIndex = '/cotton/Liuzhenping/Pan-genome/AD1-genomes/HC04/HC04.chr.fa'
bwaIndex = '/public/home/weizhang/Call_SNPS/bwa_index/HC04.chr.fa'
rule all:
	input:
		expand("{bwaIndex}.amb",bwaIndex=[bwaIndex])

#构建参考基因索引
rule bwa_index:
	input:
       		refera_fasta="{bwaIndex}"
	output:	
		"{bwaIndex}.amb"
	shell:
		""" 
		module load SAMtools/1.9
		module load BWA/0.7.17
	    bwa index {input.refera_fasta}
		"""





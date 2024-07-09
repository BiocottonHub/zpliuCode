#sampleName='S89'
# outPath='/public/home/weizhang/Call_SNPS/bams'
# inputPath='/public/home/weizhang/Call_SNPS

#S31

rule all:
	input:
		expand("/public/home/weizhang/Call_SNPS/bam_V2/{sampleName}/{sampleName}_srt_q20_redup.bam",sampleName=['S131','S276','S420']),
                expand("/public/home/weizhang/Call_SNPS/bam_V2/{sampleName}/{sampleName}_metrics.txt",sampleName=['S131','S276','S420']),
		expand("/public/home/weizhang/Call_SNPS/bam_V2/{sampleName}/{sampleName}_srt_q20_redup.bam.bai",sampleName=['S131','S276','S420'])			


#去除PCR重复
rule dedup:
	input:
		"/public/home/weizhang/Call_SNPS/bam_V2/{sampleName}/{sampleName}_sorted_q20.bam"
	output:
		bam="/public/home/weizhang/Call_SNPS/bam_V2/{sampleName}/{sampleName}_srt_q20_redup.bam",
    		metrics="/public/home/weizhang/Call_SNPS/bam_V2/{sampleName}/{sampleName}_metrics.txt",
		bai="/public/home/weizhang/Call_SNPS/bam_V2/{sampleName}/{sampleName}_srt_q20_redup.bam.bai"
	threads:1
	resources:
		mem_mb=50000
	shell:
		"""
		module load picard/2.23.9	
		java -jar /public/home/software/opt/bio/software/picard/2.23.9/picard.jar MarkDuplicates \
        	INPUT={input} \
        	OUTPUT={output.bam} \
        	METRICS_FILE={output.metrics}

		module load SAMtools/1.9
		samtools index {output.bam} {output.bai}
		"""


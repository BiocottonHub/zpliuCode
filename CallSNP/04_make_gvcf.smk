inupuSampleListFile = config["input"]
sampleList = []
with open(inupuSampleListFile, "r") as File:
    for line in File:
        sampleList.append(line.strip("\n"))


rule all:
    input:
        # vcf=expand(
        #     'gvcf/{sample}/{sample}.vcf',sample=sampleList
        # ),
        # gvcf=expand(
        #     'gvcf/{sample}/{sample}.g.vcf.gz',sample=sampleList
        # )
        # Chrvcf=expand(
        #     "mgvcf/{Chr}.g.vcf",Chr=["HC04_D12"]
        # )
        # Chr_imputedVCF=expand(
        #     "imputed/{Chr}_impute.vcf.gz",Chr=['HC04_A06',"HC04_A11"]
        # )
        # SNPs=expand(
        #     "Results/{Chr}.recode.vcf",
        #     Chr=[
        #         "HC04_A01",
        #         "HC04_A02",
        #         "HC04_A03",
        #         "HC04_A04",
        #         "HC04_A05",
        #         "HC04_A06",
        #         "HC04_A07",
        #         "HC04_A08",
        #         "HC04_A09",
        #         "HC04_A10",
        #         "HC04_A11",
        #         "HC04_A12",
        #         "HC04_A13",
        #         "HC04_D01",
        #         "HC04_D02",
        #         "HC04_D03",
        #         "HC04_D04",
        #         "HC04_D05",
        #         "HC04_D06",
        #         "HC04_D07",
        #         "HC04_D08",
        #         "HC04_D09",
        #         "HC04_D10",
        #         "HC04_D11",
        #         "HC04_D12",
        #         "HC04_D13",
        #     ],
        # ),
        # Indels=expand(
        #     "Results/{Chr}_indel.recode.vcf",
        #     Chr=[
        #         "HC04_A01",
        #         "HC04_A02",
        #         "HC04_A03",
        #         "HC04_A04",
        #         "HC04_A05",
        #         "HC04_A06",
        #         "HC04_A07",
        #         "HC04_A08",
        #         "HC04_A09",
        #         "HC04_A10",
        #         "HC04_A11",
        #         "HC04_A12",
        #         "HC04_A13",
        #         "HC04_D01",
        #         "HC04_D02",
        #         "HC04_D03",
        #         "HC04_D04",
        #         "HC04_D05",
        #         "HC04_D06",
        #         "HC04_D07",
        #         "HC04_D08",
        #         "HC04_D09",
        #         "HC04_D10",
        #         "HC04_D11",
        #         "HC04_D12",
        #         "HC04_D13",
        #     ],
        # ),
        SNPs_merge="Results/All_chr_SNP.vcf.gz",
        Indel_merge="Results/All_chr_indel.vcf.gz",


rule make_index:
    input:
        fasta="/public/home/zpliu/LZP_sQTL_A2_At/Reference_AD1/HC04.chr.fa",
    output:
        fasta_dict="/public/home/zpliu/LZP_sQTL_A2_At/Reference_AD1/HC04.chr.dict",
    resources:
        mem_mb=5000,
    threads: 1
    shell:
        """
        module load picard/2.23.9
        java -jar ${{EBROOTPICARD}}/picard.jar CreateSequenceDictionary  R={input.fasta} O={output.fasta_dict}
        """


rule bam2gvcf_in_GPU:
    # * 需要在142节点上提交任务
    input:
        bam="/public/home/weizhang/Call_SNPS/bam_V2/{sample}/{sample}_srt_q20_redup.bam",
        fasta="/public/home/zpliu/LZP_sQTL_A2_At/Reference_AD1/HC04.chr.fa",
    output:
        vcf="gvcf/{sample}/{sample}.vcf",
    resources:
        mem_mb=15000,
    threads: 1
    params:
        bindPath="/public/home/weizhang/Call_SNPSbam_V2/S96/",
    shell:
        """
        module load Singularity/3.7.3
        singularity  exec  -B {params.bindPath} --nv  $IMAGE/clara-parabricks/4.0.1-1.sif \
            pbrun haplotypecaller  --ref {input.fasta} \
            --in-bam {input.bam} \
            --out-variants {output.vcf}
        """


rule bam2gvcf_in_nonGPU:
    # * 在非GPU节点上Call SNPs
    input:
        fasta="/public/home/zpliu/LZP_sQTL_A2_At/Reference_AD1/HC04.chr.fa",
        bam="/public/home/weizhang/Call_SNPS/bam_V2/{sample}/{sample}_srt_q20_redup.bam",
    output:
        vcf="gvcf/{sample}/{sample}.g.vcf.gz",
    resources:
        mem_mb=32000,
    threads: 4
    shell:
        """
        module load GATK/4.1.9.0
        gatk --java-options '-Xmx500g -XX:ConcGCThreads={threads} -XX:ParallelGCThreads={threads}' HaplotypeCaller \
            -ERC GVCF \
            -stand-call-conf 30 \
            --native-pair-hmm-threads {threads} \
            -R {input.fasta} \
            -I {input.bam} \
            -O {output.vcf}
        """


rule merge_gvcf:
    input:
        gvcfListFile="/public/home/zpliu/LZP_sQTL_A2_At/Call_SNPs/gvcf.list",
        fasta="/public/home/zpliu/LZP_sQTL_A2_At/Reference_AD1/HC04.chr.fa",
    output:
        Chr_gvcf="mgvcf/{Chr}.g.vcf",
    threads: 20
    resources:
        mem_mb=150000,
    shell:
        """
        module load GATK/4.1.9.0
        gatk --java-options '-Xmx500g -XX:ConcGCThreads={threads} -XX:ParallelGCThreads={threads}' CombineGVCFs \
            -R  {input.fasta} \
            --variant  {input.gvcfListFile} \
            --intervals {wildcards.Chr} -O {output.Chr_gvcf}
        """


rule Call_SNPS:
    input:
        Chr_gvcf="mgvcf/HC04_D08.g.vcf",
        fasta="/public/home/zpliu/LZP_sQTL_A2_At/Reference_AD1/HC04.chr.fa",
    output:
        VCF_out="VCF/HC04_D08.vcf.gz",
    threads: 2
    resources:
        #* 20GB and run 24h
        mem_mb=20000,
    shell:
        """
        module load GATK/4.1.9.0
        gatk --java-options '-Xmx500g -XX:ConcGCThreads={threads} -XX:ParallelGCThreads={threads}' GenotypeGVCFs \
            -R {input.fasta} \
            -V  {input.Chr_gvcf}\
            --allow-old-rms-mapping-quality-annotation-data \
            -O  {output.VCF_out}
        """


rule imputed_VCF:
    input:
        rawVCF="/public/home/zpliu/LZP_sQTL_A2_At/Call_SNPs/VCF/{Chr}.vcf.gz",
    output:
        missratio_vcf=temp("imputed/{Chr}_vcftools_filter.vcf"),
        imputedVCF="imputed/{Chr}_impute.vcf.gz",
    params:
        prefix_vcf="imputed/{Chr}_impute",
    threads: 10
    resources:
        mem_mb=30000,
    shell:
        """
        module unload Java/11.0.8
        module load VCFtools/0.1.16
        module load beagle/4.1-Java-1.8.0_92
        module load BCFtools/1.8
        #* 有个别位点genotype为.,需要替换为./.
        vcftools  --gzvcf {input.rawVCF} \
            --min-alleles 2 \
            --max-alleles 2 \
            --max-missing 0.8 \
            --recode-INFO-all \
            --recode \
            -c |bcftools annotate --remove FILTER,INFO,^FORMAT/GT,FORMAT/AD - -O v \
        |sed 's/\\t\.:/\\t\.\/\.:/g' >{output.missratio_vcf}
        #* imputed model with GT 
        java -jar $EBROOTBEAGLE/beagle.jar \
            nthreads={threads} gt={output.missratio_vcf} \
            window=150 overlap=20 out={params.prefix_vcf}
        """


rule split_SNP_Indel:
    input:
        imputed_VCF="imputed/{Chr}_impute.vcf.gz",
    output:
        SNPs="Results/{Chr}.recode.vcf",
        Indels="Results/{Chr}_indel.recode.vcf",
    params:
        indel_prefix="Results/{Chr}_indel",
        snp_prefix="Results/{Chr}",
    threads: 1
    resources:
        mem_mb=10000,
    shell:
        """
        module load VCFtools/0.1.16
        vcftools --gzvcf {input.imputed_VCF} \
            --min-alleles 2 --max-alleles 2 \
            --recode --recode-INFO-all  \
            --keep-only-indels  --out {params.indel_prefix}
        vcftools --gzvcf {input.imputed_VCF} \
            --min-alleles 2 --max-alleles 2 \
            --recode --recode-INFO-all  \
            --remove-indels  --out {params.snp_prefix}
        """


rule merge_Chr_VCF:
    input:
        SNP_vcffile_list="Results/SNP_vcf.list",
        Indel_vcffile_list="Results/Indel_vcf.list",
    output:
        SNPs="Results/All_chr_SNP.vcf.gz",
        Indel="Results/All_chr_indel.vcf.gz",
    threads:2
    resources:mem_mb=20000
    shell:
        """
        module load VCFtools/0.1.16
        vcf-concat -f {input.Indel_vcffile_list} |gzip >{output.Indel}
        vcf-concat -f {input.SNP_vcffile_list} |gzip >{output.SNPs}
        """


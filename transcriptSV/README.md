## transcriptSV

Detect the SVs(indel PAVs) in the two genome feature, for example genes or transcripts.

shortcut

<img src="https://ftp.bmp.ovh/imgs/2020/12/558d13fa9b93b8c3.png" alt="show transcript sv" style="zoom:80%;" />

#### 1.dependence：

+ python3
+ [`muscle` (Multiple sequence alignment)](http://www.drive5.com/)

> :warning:	note: In the script `muscle.py`，will to call `muscle` module. what's more I have set the path be `~/software/muscle3.8.31_i86linux64`，so if  you must set the path to `muscle` to guarante script work.

#### 2.Usage:

+ `homolog` the homologous of two genome
+ `fasta1`  A genome transcripts cDNA sequence 
+ `fasta2`  B genome transcripts cDNA sequence 
+ `RNAseq1` A genome transcriptsexpression
+ `RNAseq2` B genome transcriptsexpression
+ `out`  output file

```bash
#for example
python ../transcriptSvs.py  -homolog  homologGeneList.txt -fasta1 A_transcript_cDNA.fa  -fasta2 B_transcript_cDNA.fa  -RNAseq1 A_transcripts_expression.txt  -RNAseq2 B_transcripts_expression.txt -transcriptSVs.txt 
```

#### 3.Output

In the last two columns , indicate the svs length in the two transcripts or genes.

for example:

> #0	198,5 
>
> the B transcripts has a 198bp delete and 5bp delete when compared to A transcripts

```bash
#geneid1	geneid2	isoforms1	isoform2	readCount1	readCount2	SV1	SV2
evm.TU.Ga01G0004	Ghir_A01G000040	PB.7404.1	PB.5347.1	9.0	11.0	0	4
evm.TU.Ga01G0008	Ghir_A01G000080	PB.7406.2	PB.5348.1	1.0	11.0	1,0	10,208,79,7,96,347,31,52
evm.TU.Ga01G0010	Ghir_A01G000100	PB.7407.1	PB.5351.1	9.0	8.0	0	198,5
evm.TU.Ga01G0012	Ghir_A01G000120	PB.7409.1	PB.5352.1	1.0	1.0	274,86,115	1,1
evm.TU.Ga01G0015	Ghir_A01G000150	PB.7411.3	PB.5353.1	3.0	1.0	2,749,91,0	1

```

#### 4.InputFile

+ homolog Gene List File

> only the 3th and 4th column are gene ID, other can be anythings

```bash
Aligment_98	BlockCount_711	evm.TU.Ga01G0004	Ghir_A01G000040
Aligment_98	BlockCount_711	evm.TU.Ga01G0005	Ghir_A01G000050
Aligment_98	BlockCount_711	evm.TU.Ga01G0007	Ghir_A01G000060
```

+ transcripts expression file

> + gene Id
> + transcripts Id 
> + any
> + any
> + RNA-seq  FPKM
> + PacBio full-length read count 
> + PacBio no-full-length read count 

```bash 
evm.TU.Ga01G0004	PB.7404.1	1782	2323	9.532540	9	4
evm.TU.Ga01G0004	PB.7404.2	1782	3054	0.526973	1	1
evm.TU.Ga01G0007	PB.7405.1	441	664	34.466625	1	0
evm.TU.Ga01G0008	PB.7406.2	903	2920	0.491072	1	3
```


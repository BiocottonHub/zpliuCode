### Scripts Annotion

> There are some scripts to operate GFF annotion 



#### 1.readFastaFromBed.py

paraments:

+ fasta sequence File
+ out @dict

#### 2.readgtf.py

functions:

+ `getGeneInfo`

> with a gtf input file and get a dict which contained all gene objects

@Object methods：

+ `getMaxLengthTranscript`
+ `getTranscriptNum`
+ `getTranscriptCoordinate`
+ `getTranscriptName`
+ `getChain`
+ `getAllspliceSite`
+ `getAlltranscriptsObject`

`getTranscriptInfo`

> with a gtf input file and get a dict which contained all transcript objects

methods:

+ `getLength`
+ `getgeneName`
+ `getTranscriptName`
+ `getcDNAlength`
+ `getExonCoordinate`
+ `getIntronCoordinate`
+ `getcDNAsequence`
+ `getCDSsequence`




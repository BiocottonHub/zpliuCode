import argparse
parser=argparse.ArgumentParser(description="intersect Bed with same reads")
parser.add_argument("-re",help="reference align bed")
parser.add_argument("-nore",help="noreference align bed")
parser.add_argument("-out",help="output file")
args=parser.parse_args()

reference=args.re
noreference=args.nore
out=args.out


  
###得到唯一匹配的行
def filterBed(fileName):
  with open(fileName,'r') as File1:
    lines=File1.readlines()
  print("读取文件完毕")
  tmpArry3=[]
  for line in lines:
    line=line.strip("\n").split("\t")
    if  int(line[2])-int(line[1])>=100 or int(line[2])-int(line[1])<=-100 :
      tmpArry3.append([line[0],line[1],line[2],line[3]])

  # tmpArry=[item[3] for item in lines]
  # tmpArry.sort()
  # print("获取唯一匹配的read name")
  # tmpArry2=[item for item in tmpArry if tmpArry.count(item)==1]
  # print("提取对应的行\n")
  # tmpArry4=[line for line in lines for item in tmpArry2 if line[3]==item]
  # print("进行指标过滤\n")
  seqname=list(set([line[3].split("/")[0] for line in tmpArry3]))
  print(seqname[0])
  tmpdict={}
  for line in tmpArry3:
    tmpdict[line[3]]=[line[0],line[1],line[2]]
  return [tmpdict,seqname]

### 对坐标取并集
def mergeLocal(read1,read2):
  if(read1[0]==read2[0]):
    tmp=[int(read1[1]),int(read1[2]),int(read2[1]),int(read2[2])]
    tmp.sort()
    return[read1[0],str(tmp[0]),str(tmp[-1])]
  else:
    return False

print("过滤reference bed\n")
referenceData=filterBed(reference)
print("过滤noreference bed\n")
noreferenceData=filterBed(noreference)

intesectSeqName=[seq1 for seq1 in referenceData[1] for seq2 in noreferenceData[1] if seq1==seq2] ##获取都有的序列名字

with open(out,'w') as File1:
  print("write to file\n")
  for name in intesectSeqName:
    Read1Name=name+"/1"
    Read2Name=name+"/2"
    if(Read1Name in referenceData[0] and Read2Name in referenceData[0] and Read1Name in noreferenceData[0] and Read2Name in noreferenceData[0]):
      Relocation=mergeLocal(referenceData[0][Read1Name],referenceData[0][Read2Name])
      noRelocation=mergeLocal(noreferenceData[0][Read1Name],noreferenceData[0][Read2Name])
      if(Relocation and noRelocation): ##保证坐标合并是在同一条染色体上的
        File1.write("\t".join(Relocation)+"\t"+"\t".join(noRelocation)+"\n")
      else:
        continue
    elif Read1Name not in referenceData[0] and Read2Name in referenceData[0] and Read1Name in noreferenceData[0] and Read2Name in noreferenceData[0]:
      Relocation=referenceData[0][Read2Name]
      noRelocation=mergeLocal(noreferenceData[0][Read1Name],noreferenceData[0][Read2Name])
      if(noRelocation): ##保证坐标合并是在同一条染色体上的
        File1.write("\t".join(Relocation)+"\t"+"\t".join(noRelocation)+"\n")
      else:
        continue
    elif Read1Name  in referenceData[0] and Read2Name not in referenceData[0] and Read1Name in noreferenceData[0] and Read2Name in noreferenceData[0]:
      Relocation=referenceData[0][Read1Name]
      noRelocation=mergeLocal(noreferenceData[0][Read1Name],noreferenceData[0][Read2Name])
      if(noRelocation): ##保证坐标合并是在同一条染色体上的
        File1.write("\t".join(Relocation)+"\t"+"\t".join(noRelocation)+"\n")
      else:
        continue
    elif Read1Name  not in referenceData[0] and Read2Name  in referenceData[0] and Read1Name not in noreferenceData[0] and Read2Name in noreferenceData[0]:
      Relocation=referenceData[0][Read2Name]
      noRelocation=noreferenceData[0][Read2Name]
      File1.write("\t".join(Relocation)+"\t"+"\t".join(noRelocation)+"\n")

    elif Read1Name  not in referenceData[0] and Read2Name  in referenceData[0] and Read1Name  in noreferenceData[0] and Read2Name not in noreferenceData[0]:
      Relocation=referenceData[0][Read2Name]
      noRelocation=noreferenceData[0][Read1Name]
      File1.write("\t".join(Relocation)+"\t"+"\t".join(noRelocation)+"\n")

    elif Read1Name  in referenceData[0] and Read2Name  not in referenceData[0] and Read1Name  in noreferenceData[0] and Read2Name not in noreferenceData[0]:
      Relocation=referenceData[0][Read1Name]
      noRelocation=noreferenceData[0][Read1Name]
      File1.write("\t".join(Relocation)+"\t"+"\t".join(noRelocation)+"\n")

    elif Read1Name  in referenceData[0] and Read2Name  not in referenceData[0] and Read1Name  not in noreferenceData[0] and Read2Name  in noreferenceData[0]:
      Relocation=referenceData[0][Read1Name]
      noRelocation=noreferenceData[0][Read2Name]
      File1.write("\t".join(Relocation)+"\t"+"\t".join(noRelocation)+"\n")
    elif Read1Name  in referenceData[0] and Read2Name   in referenceData[0] and Read1Name  not in noreferenceData[0] and Read2Name  in noreferenceData[0]:
      Relocation=mergeLocal(referenceData[0][Read1Name],referenceData[0][Read2Name])
      noRelocation=noreferenceData[0][Read2Name]
      if(Relocation):
        File1.write("\t".join(Relocation)+"\t"+"\t".join(noRelocation)+"\n")
      else:
        continue
    elif Read1Name  in referenceData[0] and Read2Name   in referenceData[0] and Read1Name  in noreferenceData[0] and Read2Name not in noreferenceData[0]:
      Relocation=mergeLocal(referenceData[0][Read1Name],referenceData[0][Read2Name])
      noRelocation=noreferenceData[0][Read1Name]
      if(Relocation):
        File1.write("\t".join(Relocation)+"\t"+"\t".join(noRelocation)+"\n")
      else:
        continue



  

import argparse
import re 

parse=argparse.ArgumentParser(description="finde reference location")

parse.add_argument("-A",help="sequence align file ")
parse.add_argument("-B",help="contigs reference file")
parse.add_argument("-out",help="output file")

args=parse.parse_args()

Afile=args.A 
Bfile=args.B
outfile=args.out 

contigsLen={}
print("reading contigs length")
with open(Bfile,'r') as File1:
  tmp=File1.readlines()
  for i in  range(0,len(tmp)):
    tmp[i]=tmp[i].strip("\n")
    if re.match(r'^>',tmp[i]):
      contigsLen[tmp[i].strip(">")]=len(tmp[i+1])-1
    else:
      continue
print("reading align file ")
align=[]
with open(Afile,'r') as File1:
  for line in File1.readlines():
    line=line.strip("\n").split("\t")
    align.append(line)

def getLocation(contigLoca1,contigLoca2,contigLen):
  if(contigLoca1<contigLoca2):
    return [contigLoca1-1,contigLen-contigLoca2]
  else:
    return [contigLen-contigLoca1,contigLoca2-1]

print("write out to disk")
with open(outfile,'w') as File1:
  for item in align:
    tmp=getLocation(int(item[4]),int(item[5]),contigsLen[item[3]])
    if(int(item[1])<int(item[2]) and int(item[1])-tmp[0]>=0):
      item[1]=str(int(item[1])-tmp[0])
      item[2]=str(int(item[2])+tmp[1])
      File1.write("\t".join(item)+"\n")
    elif int(item[1])<int(item[2]) and int(item[1])-tmp[0]<0:
      item[1]='0'
      item[2]=str(int(item[2])+tmp[1])
      File1.write("\t".join(item)+"\n")
    elif int(item[1])>int(item[2]) and int(item[2])-tmp[1]>=0:
      tmp2=item[2]
      item[2]=str(int(item[1])+tmp[0])
      item[1]=str(int(tmp2)-tmp[1])
      File1.write("\t".join(item)+"\n")      
    elif int(item[1])>int(item[2]) and int(item[2])-tmp[1]<0:
      tmp2=item[2]
      item[2]=str(int(item[1])+tmp[0])
      item[1]='0'
      File1.write("\t".join(item)+"\n")



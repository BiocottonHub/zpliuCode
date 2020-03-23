import argparse
parse=argparse.ArgumentParser(description="merge promter and gene body file")
parse.add_argument("-gene",help="gene bed file")
parse.add_argument("-promter",help="promter bed file")
parse.add_argument("-out",help="out  file")
args=parse.parse_args()
genefile=args.gene
promterfile=args.promter
outfile=args.out
outdict={}
with open(genefile,'r') as File1:
  genetmp=File1.readlines()
with open(promterfile,'r') as File1:
  promtertmp=File1.readlines()

for line in genetmp:
    line=line.strip("\n").split("\t")
    key="\t".join(line[0:6])
    if key not in outdict:
      outdict[key]=[[],[]]
    else:
      continue
for line in genetmp:
    line=line.strip("\n").split("\t")
    key="\t".join(line[0:6])
    outdict[key][0].append(line[6])

for line in promtertmp:
    line=line.strip("\n").split("\t")
    key="\t".join(line[0:6])
    key="\t".join(line[0:6])
    outdict[key][1].append(line[6])

with open(outfile ,'w') as File1:
  for key in outdict:
    File1.write(key+"\t"+"/".join(outdict[key][0])+"\t"+"/".join(outdict[key][1])+"\n")

import sys
import re
import argparse
'''
Usage:
  python $0 collinearitFile TADFile1 TADFile2 OutFile

共线性文件中第一部分的坐标与TAD1文件中的坐标对应
共线性文件中第二部分坐标与TAD2文件中坐标对应
'''
parser=argparse.ArgumentParser(description="Accodring genome' collinearity file to find conserve TAD Boundary in giving two genome ")
parser.add_argument("-c",help="genome collinearity file")
parser.add_argument("-TAD1",help="A genome's TAD file")
parser.add_argument("-TAD2",help="B genome's TAD file")
parser.add_argument("-out",help="output file")
parser.add_argument("-window",help="Boundary float size")
args=parser.parse_args()

collinearityFile=args.c 
TAD1=args.TAD1
TAD2=args.TAD2
outfile=args.out
window=args.window



'''
定义染色体坐标类
方法 intersectionCoords:
  判断两个坐标是否存在交集
'''


class coords:
  def __init__(self, chromose, start, end):
    self.chromose = chromose
    self.start = start
    self.end = end
    self.intersectCore=None

  def intersectionCoords(self, coordB):
    if self.chromose != coordB.chromose:
      return False
    elif self.start > coordB.end or self.end < coordB.start:
      return False
    else:
      distance=abs((self.start+self.end-coordB.start-coordB.end)/2)
      self.intersectCore=distance
      return True
  def printLocation(self):
    print(self.chromose+"\t"+str(self.start)+"\t"+str(self.end))


'''
定义共线性类
'''


class collinearity:
  def __init__(self, block1, block2):
    self.block1 = block1
    self.block2 = block2


'''
定义每个TAD的boundary区域,和滑动窗口大小
考虑左右边界,在TAD文件中右边界肯定比左边界大

'''


class TADBoundary:
  def __init__(self, TADcoords, TADName, windowSize):
    self.TADName = TADName
    self.TADcoords = TADcoords
    if(TADcoords.start==0):
      self.Boundarit = coords(TADcoords.chromose, 0, windowSize)
    else:
      self.Boundarit=coords(TADcoords.chromose,TADcoords.end-windowSize,TADcoords.end+windowSize)


'''
将共线性区段存进数组内
输入文件格式:
Chr01 4139	5266	Chr01	29416	30542
Chr01 4139	5266	Chr01	29416	30542
Chr01 4139	5266	Chr01	29416	30542
'''
collinearityArray = []
with open(collinearityFile, 'r') as collinearityFile:
  for line in collinearityFile.readlines():
    if re.match(r'#', line):
      print("异常block区域:"+line.strip("\n"))
      continue
    else:
      line = line.strip("\n").split("\t")
      if(int(line[1]) > int(line[2])):  # 染色体倒位现象
        block1 = coords(line[0], int(line[2]), int(line[1]))
      else:
        block1 = coords(line[0], int(line[1]), int(line[2]))
      if(int(line[4]) > int(line[5])):  # 染色体倒位现象
        block2 = coords(line[3], int(line[5]), int(line[4]))
      else:
        block2 = coords(line[3], int(line[4]), int(line[5]))
      collinearityArray.append(collinearity(block1, block2))

'''
将D5的TAD坐标存进boundary数组
输入文件格式:
Chr01   980000  1020000 D5_TAD_Chr01_0_1000000_1
Chr01   14130000        14170000        D5_TAD_Chr01_13650000_14150000_22
Chr01   18780000        18820000        D5_TAD_Chr01_17800000_18800000_29
'''
####初始窗口大小
initWindosSize=int(window)
TADBoundary1 = []
TADBoundary2 = []

with open(TAD1, 'r') as TADFile1:
  line1=TADFile1.readlines()
  for i in range(0,len(line1)):
    line1[i] = line1[i].strip("\n").split("\t")

with open(TAD2, 'r') as TADFile2:
  line2=TADFile2.readlines()
  for i in range(0,len(line2)):
    line2[i] = line2[i].strip("\n").split("\t")


def adjustTAD2Boundarit(line,windowSize):
  TmpTADBoundary=[]
  for i in range(0,len(line)):
    if(int(line[i][1])==0):
      TADcoords=coords(line[i][0],0,0)
      TmpTADBoundary.append(TADBoundary(TADcoords, line[i][3]+"_left", windowSize)) ##从TAD开始往共线性区域靠近
    if(i>=1 and int(line[i-1][2])!=int(line[i][1]) and int(line[i][1])!=0):
      TADcoords = coords(line[i][0], int(line[i][1]), int(line[i][1]))
      TmpTADBoundary.append(TADBoundary(TADcoords, line[i][3]+"_left", windowSize))  # 上下滑动150kb
    TADcoords = coords(line[i][0], int(line[i][2]), int(line[i][2])) ##防止0的情況
    TmpTADBoundary.append(TADBoundary(TADcoords, line[i][3], windowSize))  # 上下滑动150kb
  return TmpTADBoundary

'''
以2M的区间进行滑动，将block进行一个扫描，去除跨度比较大的block区间
'''
def deletErrorBlock(blockArray):
  deletBlock=[]
  for i in blockArray:
    length=len(blockArray)
    item=0
    for j in blockArray:
      if(j.start-i.start>1000000 or j.start-i.start <-1000000 ):
        item+=1
    if(item/length>=0.5 and item/length!=1): ##超过50%的坐标都与这个坐标不相同
       deletBlock.append(i)
  for i in deletBlock:
    blockArray.remove(i)
  if len(blockArray)!=0:
    return blockArray
  else:
    for i in deletBlock:
      print("复杂的Block区域:")
      i.printLocation()
    return False


'''
寻找150KB范围内唯一converse的boundary和不converse的bundary

'''

def findConverseBoundary(TAD1,initWindosSize):
  tmp=[]
  tmpArray=[var for var in collinearityArray if var.block1.intersectionCoords(TAD1.Boundarit)] ##获得与TAD1属于同一个基因组的collinearity类
  # 获取与TAD1距离最近的共线性区段
  if len(tmpArray)>1:
    BestCollinearity=tmpArray[0]
    for item in range(1,len(tmpArray)):
      if(tmpArray[item].block1.intersectCore<BestCollinearity.block1.intersectCore):
        BestCollinearity=tmpArray[item]
  else:
    return tmp
  tmpArry2=[var.block2 for var in tmpArray] ##获取另一个基因组对应的共线性区域，并对异常的block区域进行筛选
  filterBlock=deletErrorBlock(tmpArry2) ##获得共线性数组
  if filterBlock:
    tmpBlock=[var.start for var in filterBlock]
    tmpBlock.sort()
    start=tmpBlock[0]
    tmpBlock=[var.end for var in filterBlock]
    tmpBlock.sort()
    end=tmpBlock[-1]
    filterCoodrs=coords(filterBlock[0].chromose,start,end)
    TADBoundary2=adjustTAD2Boundarit(line2,0) ##滑动第二个TAD区间,找到尽可能多的保守TAD
    tmp=[var for var in TADBoundary2 if var.Boundarit.intersectionCoords(filterCoodrs)] ##将共线性区段与TAD2进行比较 
    if len(tmp)>1:
      BestTAD2=tmp[0] #寻找最优boundary
      for item in range(1,len(tmp)):
        if tmp[item].Boundarit.intersectionCoords(BestCollinearity.block2)<BestTAD2.Boundarit.intersectionCoords(BestCollinearity.block2):
          BestTAD2=tmp[item]
      return [BestTAD2.TADName]
    elif len(tmp) ==1:
      return [tmp[0].TADName]
    else:
      return tmp
  else:
    return tmp





outData = {}
outData["conserve"] = []
windowSize=initWindosSize
while(windowSize<=initWindosSize):
  print(">>>>>>当前TAD1滑动距离窗口大小为:"+str(windowSize))
  TADBoundary1=adjustTAD2Boundarit(line1,windowSize) ##逐渐增加TAD1的窗口大小
  for item in TADBoundary1:
    if item.TADName not in [var[0].TADName for var in outData["conserve"]]: ##当前TAD1在之前的滑动窗口下没有找到过
      tmpArray=findConverseBoundary(item,windowSize)
      if len(tmpArray)!=0: ##找到结果记录对应的TAD1名字
        print("找到结果:\t"+item.TADName+"\t"+"\t".join(tmpArray)+"\n")
        outData['conserve'].append([item,tmpArray])
  windowSize+=10000


# for item in TADBoundary1:
#   tmpArray=findConverseBoundary(item)
#   if len(tmpArray)==0:
#     outData['noconserve'].append(item.TADName)
#   else:
#     outData['conserve'].append([item.TADName,tmpArray[0]])


with open(outfile,'w') as out:
  for item in outData['conserve']:
    out.write(item[0].TADName+"\t"+"\t".join(item[1])+"\n")
  for item in TADBoundary1:
    if item.TADName not in [var[0].TADName for var in outData["conserve"]]:
      out.write(item.TADName+"\t"+"None"+"\n")


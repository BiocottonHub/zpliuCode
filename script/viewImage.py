'''
Descripttion: 
version: 
Author: zpliu
Date: 2020-12-25 10:55:40
LastEditors: zpliu
LastEditTime: 2020-12-25 14:31:35
@param: 
'''
import os
import psutil  # to kill the picture process
from PIL import Image
import re
import sys

# @async1


def openPicture(fileName):
    '''
    @Descripttion:  transfer picture 
    @param: fileName  picture file name@str
    @return: 
    '''
    img = Image.open(fileName)
    img.show()
    img.close()


def setChromsomes(geneId, fileName):
    '''
    @Descripttion: set abnormal chromsomes id
    @param: geneId @str
    @param: fileName @str
    @return: 
    '''
    chromsomes = input("染色体编号:")
    chromsomes = re.split(r'\s+', chromsomes.strip("\n"))
    with open(fileName, 'a') as File:
        File.write(geneId+"\t"+",".join(chromsomes)+"\n")

    for proc in psutil.process_iter():
        if proc.name() == "Microsoft.Photos.exe":  # turn off the photo
            proc.kill()  # 关闭该proces
    return


def needViewImage(dirname, fishedFile):
    '''
    @Descripttion:  filter fished gene id 
    @param: picture dirname @str
    @param: out filename @str
    @return: 
    '''
    all = []
    tmp = []
    with open(fishedFile, 'r') as File:
        for line in File:
            line = line.strip("\n").split("\n")
            tmp.append(line[0])
    for fileName in os.listdir(dirname):
        if re.search(r'png$', fileName):
            if re.sub(r'_manhattan.*', '', fileName) in tmp:
                pass
            else:
                all.append(re.sub(r'_manhattan.*', '', fileName))
    return all


outFile = sys.argv[2]
FilePath = os.path.abspath(sys.argv[1])
outFilePath = os.path.abspath(sys.argv[2])
Allgenes = needViewImage(FilePath, outFilePath)
count = 0
for geneId in Allgenes:
    count += 1
    openPicture(os.path.join(FilePath, geneId+'_manhattan.png'))
    print("当前进度{}/{}".format(str(count), str(len(Allgenes))), end="\t")
    setChromsomes(geneId, outFilePath)
#[i.join() for i in Task]

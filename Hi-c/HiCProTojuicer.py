'''/*
 * @Author: mikey.zhaopeng 
 * @Date: 2020-07-05 10:44:49 
 * @Last Modified by:   mikey.zhaopeng 
 * @Last Modified time: 2020-07-05 10:44:49 
 */
'''


import sys
from multiprocessing import Process, Queue
import re
import gc
import os
import fcntl


class PartitionFile(object):
    def __init__(self, fileName, jobsNum):
        self.fileName = fileName
        self.blockNum = jobsNum

    def partion(self):
        fd = open(self.fileName, 'r')
        fd.seek(0, 2)  # 移动文件指针到文件尾,用于获取文件大小
        fileSize = fd.tell()  # 获取文件字符数
        Pos_list = []  # 指针坐标，数组
        blockSize = int(fileSize/self.blockNum)
        start_Pos = 0  # 文件初始指针
        for i in range(self.blockNum):
            if i == self.blockNum-1:
                end_Pos = fileSize-1  # 最后一个文件区块为文件结尾
                Pos_list.append((start_Pos, end_Pos))
                break
            end_Pos = start_Pos+blockSize-1  # 均匀分配每个区块
            Pos_list.append((start_Pos, end_Pos))
            start_Pos = end_Pos+1  # 下一个区块的开始坐标
        fd.close()
        return Pos_list


'''
创建临时文件夹
'''


def mkdir(path):
    path = path.strip()
    path = path.rstrip("\\")
    isExists = os.path.exists(path)
    if not isExists:
        os.makedirs(path)
        return True
    else:
        return False


class readProcess(Process):  # 这个括号表示继承threading.Thread类
    def __init__(self, job_name, fileName, queue, start_Pos, end_Pos, processFunction):
        super(readProcess, self).__init__()
        self.name = job_name
        self.start_Pos = start_Pos
        self.end_Pos = end_Pos
        self.fileName = fileName
        self.processFunction = processFunction
        self.queue = queue  # 用于存储染色体编号
        self.outData = {}

    def run(self):
        print("Process: " + self.name+": reading file...")
        self.reader()
        print("Process: " + self.name+": begin to write temporary data to file...")
        for key in self.outData.keys():
            self.queue.put(key, block=True, timeout=None)
            with open(path+"/"+key, 'a') as File:
                fcntl.flock(File.fileno(), fcntl.LOCK_EX)
                for item in self.outData[key]:
                    File.write(item[0]+" "+key+" "+" ".join(item[1:])+"\n")
        print("Process: " + self.name+": completed write...")

    def reader(self):
        fd = open(self.fileName, 'r')
        if self.start_Pos != 0:
            fd.seek(self.start_Pos-1)
            if fd.read(1) != '\n':  # 当前初始位置不是行首,移动到下一行行首
                fd.readline()
                self.start_Pos = fd.tell()
        fd.seek(self.start_Pos)  # 将文件指针定位到区块的行首
        while self.start_Pos <= self.end_Pos:  # 开始按行读取文件并且进行操作
            tmp = fd.readline()
            self.processFunction(tmp, self.outData)  # 将每行结果存进字典内
            self.start_Pos = fd.tell()  # 读完一行后，自动调整开始位置
        fd.close()
        return


class sortProcess(Process):
    def __init__(self, job_name, fileName):
        super(sortProcess, self).__init__()
        self.name = job_name
        self.fileName = fileName
        self.outData = {}

    def run(self):
        print("sorting chrosome: "+self.name+"...")
        self.sortFile()
        print("chrosome: "+self.name+"  ok...")

    def sortFile(self):
        with open(path+"/"+self.fileName, 'r') as File:
            for line in File.readlines():
                line = line.split(" ")
                try:
                    # [tmpDict[line[3]], line[1],line[2], '0', tmpDict[line[6]], line[4], line[5], '0']
                    # 都是同一条染色体对应的Chr1-Chr2 Chr1-Chr3
                    self.outData[line[5]].append(
                        [line[0], line[2], line[3], line[4], line[6], line[7]])
                except KeyError:
                    self.outData[line[5]] = [
                        [line[0], line[2], line[3], line[4], line[6], line[7]]]
        with open(path+"/"+self.fileName+"_sorted", 'w') as File:
            sortKey = sorted(self.outData)
            for key in sortKey:
                for item in self.outData[key]:
                    File.write(item[0]+" "+self.fileName+" " +
                               " ".join(item[1:4])+" "+key+" "+item[-2]+" "+item[-1])


def HicProTojuicer(fileName, ProcessNum, scaffoldPrefix, OutFile):
    def processFunction(line, ouput):  # 对每行文件进行处理，并将结果存进二维数组
        line = line.split("\t")
        tmpDict = {'+': '0', '-': '1'}
        if re.match(scaffoldPrefix, line[1]) or re.match(scaffoldPrefix, line[4]):
            return False
        else:
            # 每个进程处理后的数据
            try:
                # 键值不再存了，减少内存消耗
                # [tmpDict[line[3]], line[1],line[2], '0', tmpDict[line[6]], line[4], line[5], '0']
                ouput[line[1]].append([tmpDict[line[3]],
                                       line[2], '0', tmpDict[line[6]], line[4], line[5], '0'])
            except KeyError:
                ouput[line[1]] = []
                ouput[line[1]] = [[tmpDict[line[3]],
                                   line[2], '0', tmpDict[line[6]], line[4],  line[5], '0']]
    global path
    path = 'tmp'+str(int(time.time()))
    mkdir(path)
    workQueue = Queue()  # 用于存放子进程文件数据
    read_jobs = []
    sort_jobs = []
    chrosomes = []
    pos_list = PartitionFile(fileName, ProcessNum).partion()  # 存放所有文件指针坐标
    for i in range(ProcessNum):
        position = pos_list[i]
        myprocess = readProcess(
            str(i), fileName, workQueue, position[0], position[1], processFunction)
        myprocess.start()
        read_jobs.append(myprocess)
    for i in read_jobs:
        i.join()
    while True:
        try:
            chrosomes.append(workQueue.get(block=True, timeout=1))  # 获取子进程数据
        except:
            break
    for i in list(set(chrosomes)):
        myprocess = sortProcess(str(i), i)  # 排序进程
        myprocess.start()
        sort_jobs.append(myprocess)
    for i in sort_jobs:
        i.join()
    print("merge chrosomes to a single file...")
    os.system('cat '+path+"/*_sorted  >"+OutFile)
    print("completed!")
    print("there are some temporary file in directory: <./" +
          path+">\n you can remove it by yourself!")
    return


if __name__ == "__main__":
    import argparse
    import time
    parser = argparse.ArgumentParser(
        description="accroding gtf File and AS event find flanking Exon sequence Tag")
    parser.add_argument("-i", help="HiC-Pro raw data file")
    parser.add_argument("-t", help=" threads number")
    parser.add_argument("-s", help="the prefix of scaffold")
    parser.add_argument("-o", help="output file")
    args = parser.parse_args()
    start_time = time.time()
    HicProTojuicer(args.i, int(args.t), args.s, args.o)  # 主程序调试
    end_time = time.time()
    print('Cost Time is {:.2f}'.format(end_time-start_time))

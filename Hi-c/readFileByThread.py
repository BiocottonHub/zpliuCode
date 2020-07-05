'''/*
 * @Author: mikey.zhaopeng 
 * @Date: 2020-07-05 10:45:07 
 * @Last Modified by:   mikey.zhaopeng 
 * @Last Modified time: 2020-07-05 10:45:07 
 */
'''

import time
from threading import Thread
from queue import Queue
import re
import sys


class PartitionFile(object):
    def __init__(self, fileName, threadNum):
        self.fileName = fileName
        self.blockNum = threadNum

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
            # if end_Pos >= fileSize:
            #   end_Pos=fileSize-1
            # if start_Pos >= fileSize:
            #   break
            Pos_list.append((start_Pos, end_Pos))
            start_Pos = end_Pos+1  # 下一个区块的开始坐标
        fd.close()
        return Pos_list


class readThread(Thread):  # 这个括号表示继承threading.Thread类
    def __init__(self, thread_name, thread_queue, fileName, start_Pos, end_Pos, processFunction):
        super(readThread, self).__init__()
        self.name = thread_name
        self.queue = thread_queue
        self.start_Pos = start_Pos
        self.end_Pos = end_Pos
        self.fileName = fileName
        self.processFunction = processFunction

    def run(self):
        print("线程" + self.name+": 开始读取文件...")
        self.reader()
        print("线程" + self.name+": 读取完成...")

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
            if self.processFunction(tmp):
                self.queue.put(tmp, block=True,
                               timeout=5)  # 队列满了，就阻塞进行等待5秒，还满就退出
            self.start_Pos = fd.tell()  # 读完一行后，自动调整开始位置
        fd.close()
        return


class getThread(Thread):
    def __init__(self, name, queue):
        super(getThread, self).__init__()
        self.name = name
        self.queue = queue
        self.processFunction = processFunction

    def run(self):
        while True:
            try:
                # 队列不会为空，不等待
                line = self.queue.get(block=True, timeout=1)
                out.append(line)
            except:
                # print(e)
                break

    def getData(self):
        print(self.out)
        return self.out


def defaultProcessFunction(line):  # 对行数据不做处理的默认函数
    return line


def readFileByThread(fileName, ThreadNum, out, processFunction=defaultProcessFunction):
    # 设置队列长度
    workQueue = Queue()
    # 线程池
    readThreads = []
    # getThreads = []
    pos_list = PartitionFile(fileName, ThreadNum).partion()
    for i in range(len(pos_list)):
        postion = pos_list[i]
        mythread = readThread(str(i), workQueue, fileName,
                              postion[0], postion[1], processFunction)  # 初始化线程,设置预处理函数
        mythread.start()  # 启动线程
        # getdataThread = getThread(str(i), workQueue)
        # getdataThread.start()
        readThreads.append(mythread)  # 添加到线程池
        # getThreads.append(getdataThread)  # 添加到线程池
    for i in readThreads:
        i.join()  # 等待所有线程完成
    # for i in getThreads:
    #     i.join()  # 等待所有线程完成
    #     out += i.getData()
    # out.append(workQueue.get(block=True, timeout=1))
    return


if __name__ == "__main__":
    def processFunction(line):  # 对每行文件进行处理，并将结果存进二维数组
        line = line.split("\t")
        tmpDict = {'+': '0', '-': '1'}
        if re.match('^S', line[1]) or re.match('^S', line[4]):
            return False
        else:
            return line
    out = []
    start_time = time.time()
    readFileByThread(sys.argv[1], int(sys.argv[2]), out,  processFunction)
    end_time = time.time()
    print('Cost Time is {:.2f}'.format(end_time-start_time))

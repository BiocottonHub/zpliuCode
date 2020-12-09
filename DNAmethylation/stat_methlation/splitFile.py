import time
import os
import sys
import gzip
from extract_read_count import *


def mkdir(pathStr):
    '''
    @Descripttion: make dir by name
    @param: path @str
    @return:  bool
    '''
    folder = os.path.exists(pathStr)
    if not folder:
        os.makedirs(pathStr)
        print(pathStr)
        return 0
    else:
        print(pathStr+":is exist!")
        return 1


def splitFileByCout(inputFilename, splitLines, FilePath, ):
    '''
    @Descripttion: 
    @param: 
    @return: 
    '''
    try:
        with gzip.open(inputFilename, 'r') as File:
            File.readline()  # head
            tmp = []
            fileIndex = 1
            for line in File:
                line = line.decode()
                tmp.append(line)
                if len(tmp) == splitLines:
                    with open(os.path.join(FilePath, str(fileIndex)), 'w') as File2:
                        File2.writelines(tmp)
                    tmp = []
                    fileIndex += 1
            if len(tmp) != 0:
                with open(os.path.join(FilePath, str(fileIndex)), 'w') as File2:
                    File2.writelines(tmp)
    except OSError:
        # not gizp file
        with open(inputFilename, 'r') as File:
            File.readline()  # head
            tmp = []
            fileIndex = 1
            for line in File:
                tmp.append(line)
                if len(tmp) == splitLines:
                    with open(os.path.join(FilePath, str(fileIndex)), 'w') as File2:
                        File2.writelines(tmp)
                    tmp = []
                    fileIndex += 1
            if len(tmp) != 0:
                with open(os.path.join(FilePath, str(fileIndex)), 'w') as File2:
                    File2.writelines(tmp)


if __name__ == "__main__":
    tmp = sys.argv[2]+str(int(time.time()))
    if mkdir(tmp):
        exit()
    splitFileByCout(sys.argv[1], 100000000, tmp)
    sys.exit(tmp)

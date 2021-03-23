'''
Descripttion: 
version: 
Author: zpliu
Date: 2021-03-23 20:13:29
LastEditors: zpliu
LastEditTime: 2021-03-23 21:01:05
@param: 
'''
'''
@usage:

python $0 BAM_file feature_file.txt outFile.txt process number
'''




from utils.get_region_PairedReadCount import getReadcount
from multiprocessing import Pool
import pandas as pd
import numpy as np
import sys
def assagmentRegion(inputData: object):
    chromsome, start, end, stand, BamFile, featureId = inputData
    return getReadcount(chromsome, start, end, stand, BamFile, featureId)


def getAllfeature(featureFile: str):
    '''
    @@Descripttion: get 
    @@param: chr start end stand featureId
    @@return: pandas.core.dataFrame
    '''
    dataFrame = pd.read_csv(featureFile, sep="\t", header=0, index_col=None)
    return dataFrame


def calculateReadCount(BamFile, FeatureFile, processNum, outFile):
    dataframe = getAllfeature(FeatureFile)
    p = Pool(processNum)
    poolJobs = []
    # * i : index,rowdata
    for splitList in np.array_split(dataframe.values, processNum):
        splitList = [np.insert(i, -1, BamFile) for i in splitList]
        res = p.map(assagmentRegion, splitList)
        poolJobs.append(res)
    # * close the pool
    p.close()
    p.join()
    # print(poolJobs[5])
    # * merge all data from multiprocess
    result = np.vstack([poolJobs[i] for i in range(processNum)])
    with open(outFile, 'w') as File:
        File.write("chromsome\tstart\tend\tstand\tfeatureId\treadCount\n")
        for item in result:
            File.write("\t".join([str(i) for i in item])+"\n")


if __name__ == "__main__":
    BamFile = sys.argv[1]
    FeatureFile = sys.argv[2]
    outFile = sys.argv[3]
    processNum = int(sys.argv[4])
    calculateReadCount(BamFile, FeatureFile, processNum, outFile)

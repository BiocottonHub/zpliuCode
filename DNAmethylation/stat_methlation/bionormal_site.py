'''
Descripttion:
version:
Author: zpliu
Date: 2020-12-07 10:11:08
LastEditors: zpliu
LastEditTime: 2020-12-10 12:49:30
@param:python bionormal.py CpG_read.gz CpG_bionormalOutFile.txt processNum
'''
import multiprocessing
from scipy import stats
import argparse


def BioTest(AllData, processId):
    '''
    @Descripttion:
    @param: AllData: all Data site @dict
    @param TestKey: all test line
    @return:
    '''
    output = ()
    for item in AllData:
        item = item.strip("\n").split("\t")
        count = int(item[2])+int(item[3])
        p_value = stats.binom.pmf(int(item[2]), count, 0.006)
        output += ("\t".join(item)+"\t"+str(p_value)+"\n",)  # this is a tuple
    return output


def accumulateResult(callbackData):
    '''
    @Descripttion:
    @param:
    @return:
    '''
    global outArray
    outArray += callbackData
    return None


def processBionormal(AllSiteLine: dict, processNum: int):
    '''
    @Descripttion: multprocesses calculate Bionormal
    @param: AllSiteLines @dict
          {'Chr12*local'@str:[souportReadCount@int,nosouportReadcount@int]}
    @return: None
    '''
    process_average_line = int(len(AllSiteLine)/processNum)
    p = multiprocessing.Pool(processNum)
    for processId in range(0, processNum):
        if processId == processNum-1:
            start = processId*process_average_line
            end = len(AllSiteLine)
        else:
            start = processId*process_average_line
            end = (processId+1)*process_average_line
        # multprocess Bionormal
        p.apply_async(BioTest,
                      (AllSiteLine[start:end], processId),
                      callback=accumulateResult)
    p.close()
    p.join()
    return None


if __name__ == "__main__":
    # 创建参数对象
    parser = argparse.ArgumentParser(
        description="extract Bionormal ")
    parser.add_argument("-process", help="process num")
    parser.add_argument("-input", help="input file")
    parser.add_argument("-output", help="output file")
    args = parser.parse_args()
    with open(args.input, 'r') as File:
        AllSiteLine = File.readlines()
    outArray = []
    processBionormal(AllSiteLine, int(args.process))
    with open(args.output, 'w') as File:
        for item in outArray:
            File.write(item)

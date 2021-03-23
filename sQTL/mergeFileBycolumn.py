'''
Descripttion:
version:
Author: zpliu
Date: 2021-03-10 23:07:24
LastEditors: zpliu
LastEditTime: 2021-03-11 11:03:39
#! @param: str@merge file
#! @param: str@add file
#! @param: str@column name of add file
#! @param: str@column id
#! @param: str@column id
'''
import sys
import pandas as pd
import os


def readTable(Filename, columnList: list, headerName: str):

    # ? Add a new file and extract the columns
    # *
    data = pd.read_csv(Filename, sep="\t",
                       usecols=columnList, header=0)
    #! reset the last columen name
    retList = list(data.columns)
    # the column will be merge
    retList[-1] = headerName
    data.columns = retList
    return data


def mergeTable(outFileName: str, addFile: str, readcols: list, headerName: str):
    '''
    # * merge Table file 
    # *
    # *
    '''
    try:
        fixtFileSize = os.path.getsize(outFileName)
    except FileNotFoundError:
        fixtFileSize = None
    if fixtFileSize:
        # fixedTable = readTable(fixedTable, [], [])
        # ? get the pickle file
        fixedTable = pd.read_pickle(outFileName)
        #! only get the  information
        addFileTable = readTable(addFile, readcols, headerName)[headerName]
        #! merge the fixed file with  the last column
        # ? merge according location and information will get None if without the row
        mergeTable = pd.concat([fixedTable, addFileTable], axis=1)
        # mergeTable.to_pickle(outFileName, index=False, header=True)
        mergeTable.to_pickle(outFileName)
        return
    else:
        # ? get all columns when fixed file is empty
        addFileTable = readTable(addFile, readcols, headerName)
        # * pickle the outfile
        addFileTable.to_pickle(outFileName)
        return


if __name__ == "__main__":
    fixedFile = sys.argv[1]
    addFile = sys.argv[2]
    headerName = sys.argv[3]
    columnIdx = [int(i) for i in sys.argv[4:]]
    mergeTable(fixedFile, addFile, columnIdx, headerName)

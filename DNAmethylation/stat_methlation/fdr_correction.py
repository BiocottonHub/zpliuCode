'''
Descripttion: 
version: 
Author: zpliu
Date: 2020-12-12 11:37:47
LastEditors: zpliu
LastEditTime: 2020-12-12 21:21:20
@param: 
'''
import sys
import pandas as pd
import numpy as np
import statsmodels.stats.multitest as multi


def fdr_coorect(inputFile, outFile):
    '''
    @Descripttion:  p-value correct 
    @param: 
    @return: 
    '''
    data = pd.read_csv(inputFile, sep="\t", header=None, low_memory=False)
    # drop some lines with SRR start
    dropRow = data[data[0].str.startswith('SRR')]
    data.drop(dropRow.index, inplace=True)
    # mult corrected p-value and add column to data.frame
    data.loc[:, 'fdr'] = multi.multipletests(np.array(data[4]))[1]
    # add the same location to data.frame
    data.loc[:, 5] = data[1]
    # adjust the order of data.frame
    data = data[[0, 1, 5, 2, 3, 4, 'fdr']]
    data.to_csv(outFile, sep="\t", header=None, index=None)
    return None


if __name__ == "__main__":
    fdr_coorect(sys.argv[1], sys.argv[2])

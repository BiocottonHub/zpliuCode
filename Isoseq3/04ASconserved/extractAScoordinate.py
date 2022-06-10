'''
Descripttion: 
version: 
Author: zpliu
Date: 2021-01-04 10:00:14
LastEditors: zpliu
LastEditTime: 2021-02-16 08:24:37
@param: 
'''

import re
import sys


def getEventLength(ASevent):
    '''
    @Descripttion:  according AS file to get AS coordinate
    @param: ASevent @str
    @return: AS message @dict
    '''
    if re.search('RI', ASevent):
        geneid = ASevent.split(";")[0]
        local = re.search(r':([0-9]*)-([0-9]*):', ASevent)
        return {
            'length': abs(int(local.group(1))-int(local.group(2))),
            'coordinate': [int(local.group(1)), int(local.group(2))-1],
            'event':  geneid+";"+'RI'+":"
            + ASevent.split(":")[1]+":"+str(int(local.group(1)))+"-" +
            str(int(local.group(2))-1)+":"+ASevent[-1],
            'chromsome': ASevent.split(":")[1],
            'stand': ASevent[-1]
        }
    elif re.search('SE', ASevent):
        local = re.search(r'-([0-9]*):([0-9]*)-', ASevent)
        return {
            'length': abs(int(local.group(1))-int(local.group(2))),
            'coordinate': [int(local.group(1))-1, int(local.group(2))],
            'event': ASevent,
            'chromsome': ASevent.split(":")[1],
            'stand': ASevent[-1]
        }
    elif re.search('A3', ASevent):
        local = re.search(r':([0-9]*)-([0-9]*):([0-9]*)-([0-9]*):', ASevent)
        tmp = sorted([int(local.group(1)), int(local.group(2)),
                      int(local.group(3)), int(local.group(4))])
        return {
            'length': tmp[-1]-tmp[0],
            'coordinate': [tmp[0], tmp[-1]],
            'event': ASevent,
            'chromsome': ASevent.split(":")[1],
            'stand': ASevent[-1]
        }
    elif re.search('A5', ASevent):
        local = re.search(r':([0-9]*)-([0-9]*):([0-9]*)-([0-9]*):', ASevent)
        tmp = sorted([int(local.group(1)), int(local.group(2)),
                      int(local.group(3)), int(local.group(4))])
        return {
            'length': tmp[-1]-tmp[0],
            'coordinate': [tmp[0], tmp[-1]],
            'event': ASevent,
            'chromsome': ASevent.split(":")[1],
            'stand': ASevent[-1]
        }
    else:
        local = re.search(r':([0-9]*)-([0-9]*)', ASevent)
        tmp = [int(local.group(1)), int(local.group(2))]
        geneid = ASevent.split(";")[0]
        return {
            'length': tmp[-1]-tmp[0],
            'coordinate': [tmp[0], tmp[-1]],
            'event': geneid+";"+'CI'+":"
            + ASevent.split(":")[1]+":"+str(tmp[0])+"-" +
            str(tmp[1])+":"+ASevent[-1],
            # 'chromsome': ASevent.split(":")[0].split(";")[-1]
            'chromsome': ASevent.split(":")[1],
            'stand': ASevent[-1]
        }


if __name__ == "__main__":
    with open(sys.argv[3], 'w') as File:
        with open(sys.argv[1], 'r') as File1:
            for line in File1:
                line = line.strip("\n").split("\t")
                out = getEventLength(line[int(sys.argv[2])-1])
                File.write(out['chromsome']+"\t" +
                           "\t".join([str(i)for i in out['coordinate']]
                                     )+"\t"+out['event']
                           + "\t"+"1\t"+out["stand"]+"\n")

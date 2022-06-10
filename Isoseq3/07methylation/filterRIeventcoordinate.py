'''
Descripttion: 
version: 
Author: zpliu
Date: 2021-01-04 10:00:14
LastEditors: zpliu
LastEditTime: 2021-01-11 20:46:19
@param: 
'''

import re
import sys


def filterItem(items, outfile, PRIDict):
    """
    docstring
    """
    out = {}
    for line in items:
        out[line[0]] = out.get(
            line[0], [line[1], int(line[2]), int(line[3])])
        if abs(int(line[2])-int(line[3])) < abs(out[line[0]][1]-out[line[0]][2]):
            out[line[0]] = [line[1], int(line[2]), int(line[3])]
        else:
            pass
    out2 = {}
    for key, value in out.items():
        if value[0] not in out2:
            out2[value[0]] = out2.get(value[0], [key, value[1], value[2]])
        else:
            if abs(value[1]-value[2]) < abs(out2[value[0]][1]-out2[value[0]][2]):
                out2[value[0]] = [key, value[1], value[2]]
    with open(outfile, 'w') as File:
        for key, value in out2.items():
            File.write(
                value[0]+"\t"+key+"\t" +
                str(value[1]) +
                "\t"+str(value[2])+"\t" +
                str(PRIDict[value[0]])+"\t" +
                str(PRIDict[key]) + "\n"
            )


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
            'length': abs(int(local.group(1))-int(local.group(2)))-1,
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
            'length': tmp[-1]-tmp[0]+1,
            'coordinate': [tmp[0]-1, tmp[-1]],
            'event': geneid+";"+'CI'+":"
            + ASevent.split(":")[1]+":"+str(tmp[0]-1)+"-" +
            str(tmp[1])+":"+ASevent[-1],
            # 'chromsome': ASevent.split(":")[0].split(";")[-1]
            'chromsome': ASevent.split(":")[1],
            'stand': ASevent[-1]
        }


def getCoordinate(inputfile, column, outFile):
    """
    docstring
    """
    with open(outFile, 'w') as File:
        with open(inputfile, 'r') as File1:
            for line in File1:
                line = line.split("\t")
                out = getEventLength(line[int(column)-1])
                File.write(out['chromsome']+"\t" +
                           "\t".join([str(i)for i in out['coordinate']]
                                     )+"\t"+out['event']
                           + "\t"+"1\t"+out["stand"]+"\n")


def getPRI(ASfile):
    out = {}
    with open(ASfile, 'r') as File:
        for line in File:
            line = line.strip("\n").split("\t")
            if re.search('constitutivein', line[1]):
                eventObject = getEventLength(line[1])
                out[eventObject['event']] = 0
            else:
                eventObject = getEventLength(line[1])
                out[eventObject['event']] = out.get(
                    eventObject['event'], float(line[6]))
            eventObject = getEventLength(line[0])
            out[eventObject['event']] = out.get(
                eventObject['event'], float(line[4]))
    return out


def filterEvent(inputfile,  outFile):
    out = []
    with open(inputfile, 'r') as File1:
        for line in File1:
            line = line.split("\t")
            out1 = getEventLength(line[0])
            out2 = getEventLength(line[1])
            out.append([out1['event'],
                        out2['event'],
                        out1['length'],
                        out2['length']
                        ])
    PRIDict = getPRI(inputfile)

    filterItem(out, outFile, PRIDict)


if __name__ == "__main__":
    filterEvent(sys.argv[1], sys.argv[2])

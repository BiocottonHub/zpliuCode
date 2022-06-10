'''
Descripttion:
version:
Author: zpliu
Date: 2021-01-11 15:21:36
LastEditors: zpliu
LastEditTime: 2021-01-11 15:40:49
@param:
'''
import os
import re
import sys


def getEventLoacal(ASevent):
    '''
    @Descripttion:  according AS file to get AS coordinate
    @param: ASevent @str
    @return: AS message @dict
    '''
    if re.search('RI', ASevent):
        local = re.search(r':([0-9]*)-([0-9]*):', ASevent)
        return {
            'length': abs(int(local.group(1))-int(local.group(2)))-1,
            'coordinate': [int(local.group(1))+1, int(local.group(2))-1],
            'event': ASevent,
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
        return {
            'length': tmp[-1]-tmp[0]+1,
            'coordinate': [tmp[0], tmp[-1]+1],
            'event': ASevent,
            # 'chromsome': ASevent.split(":")[0].split(";")[-1]
            'chromsome': ASevent.split(":")[1],
            'stand': ASevent[-1]
        }


def writeFile(outFile, outItems):
    with open(outFile, 'a+') as File:
        for item in outItems:
            File.write(item)


def readeventSequence(fastaFile):
    out = {}
    with open(fastaFile, 'r') as File:
        for line in File:
            if re.match(r'^>', line):
                geneId = re.split("::", line)[0].strip(">")
            else:
                out[geneId] = ">"+geneId+"\n"+line
    return out


def muscle(infasta, softwarePath):
    out = os.popen(
        "printf \"%s\"|%s  -clw  2>/dev/null" % (infasta, softwarePath)).read()
    return len(re.findall(r'\*', out))


if __name__ == "__main__":
    eventFile = sys.argv[1]
    fastaFile = sys.argv[2]
    musclePath = sys.argv[3]
    with open(sys.argv[4], 'w') as File:
        pass
    fastaDict = readeventSequence(fastaFile)
    with open(eventFile, 'r') as File:
        for line in File:
            line = line.strip("\n").split("\t")
            eventObject1 = getEventLoacal(line[0])
            eventObject2 = getEventLoacal(line[1])
            fastasequence1 = fastaDict[line[0]]
            fastasequence2 = fastaDict[line[1]]
            infasta = fastasequence1+fastasequence2
            idnetityLength = muscle(infasta, musclePath)
            out = "\t".join((
                eventObject1['chromsome'],
                str(eventObject1['coordinate'][0]),
                str(eventObject1['coordinate'][1]),
                eventObject2['chromsome'],
                str(eventObject2['coordinate'][0]),
                str(eventObject2['coordinate'][1]),
                eventObject1['event'],
                eventObject2['event'],
                str(idnetityLength),
                str(eventObject1['length']),
                str(eventObject2['length']),
            )) + "\n"
            writeFile(sys.argv[4], out)

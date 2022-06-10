'''
Descripttion:
version:
Author: zpliu
Date: 2021-01-09 13:57:16
LastEditors: zpliu
LastEditTime: 2021-01-09 21:58:20
@param:
'''
import sys
import re
import os


def getConservedCytosins(muscleOutStr: str,
                         Absolutelocation1: tuple,
                         Absolutelocation2: tuple,
                         event1: str,
                         event2: str,
                         stand1: 'str',
                         stand2: 'str'

                         ):
    '''
    @Descripttion:
    @param:
    @return:
    '''
    out = []
    '''
    x: stand,
    y: absolutelocal
    z: slide base
    '''
    def f(x, y, z): return (y[0]+z, y[0] +
                            z) if x == "+" else (y[-1]-z, y[-1]-z)
    for flagment in muscleOutStr.split("\n\n")[1::]:
        sequenceflagment1, sequenceflagment2 = [
            re.split("\s+", i)[-1] for i in flagment.strip("\n").split("\n")[0:2]]
        '''
        score line  index with 38 gap length  string
        '''
        for sequence1, sequence2 in zip(sequenceflagment1, sequenceflagment2):
            if re.search("^-[ATCGatcgNn]", sequence1+sequence2):
                Absolutelocation1 = f(
                    stand1, Absolutelocation1, 0)  # delet a base
                Absolutelocation2 = f(stand2, Absolutelocation2, 1)
                out.append(
                    "\t".join((sequence1,
                               sequence2,
                               'None',
                               str(Absolutelocation2[0])
                               )
                              )+"\t"+event1+"\t"+event2+"\n")
            elif re.search("^[ATCGatcgNn]-", sequence1+sequence2):
                Absolutelocation1 = f(stand1, Absolutelocation1, 1)
                Absolutelocation2 = f(
                    stand2, Absolutelocation2, 0)  # delet a base
                # print(Absolutelocation2)
                out.append(
                    "\t".join((sequence1,
                               sequence2,
                               str(Absolutelocation1[0]),
                               'None'
                               )
                              )+"\t"+event1+"\t"+event2+"\n")
            else:
                Absolutelocation1 = f(stand1, Absolutelocation1, 1)
                Absolutelocation2 = f(
                    stand2, Absolutelocation2, 1)
                out.append(
                    "\t".join((sequence1,
                               sequence2,
                               str(Absolutelocation1[0]),
                               str(Absolutelocation2[0])
                               )
                              )+"\t"+event1+"\t"+event2+"\n")
    return out


def readeventSequence(fastaFile):
    out = {}
    with open(fastaFile, 'r') as File:
        for line in File:
            if re.match(r'^>', line):
                geneId = re.split("::", line)[0].strip(">")
            else:
                out[geneId] = ">"+geneId+"\n"+line
    return out


def getEventLoacal(ASevent):
    '''
    @Descripttion:  according AS file to get AS coordinate
    @param: ASevent @str
    @return: AS message @dict
    '''
    if re.search('RI', ASevent):
        local = re.search(r':([0-9]*)-([0-9]*):', ASevent)
        return {
            'length': abs(int(local.group(1))-int(local.group(2))),
            'coordinate': [int(local.group(1)), int(local.group(2))],
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
            'length': tmp[-1]-tmp[0],
            'coordinate': [tmp[0]-1, tmp[-1]+1],
            'event': ASevent,
            # 'chromsome': ASevent.split(":")[0].split(";")[-1]
            'chromsome': ASevent.split(":")[1],
            'stand': ASevent[-1]
        }


def writeFile(outFile, outItems):
    with open(outFile, 'a+') as File:
        for item in outItems:
            File.write(item)


if __name__ == "__main__":
    eventFile = sys.argv[1]
    fastaFile = sys.argv[2]
    musclePath = sys.argv[3]
    '''
    empty the out file
    '''
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
            muscleOut = os.popen(
                "printf \"%s\"|%s  -clw  2>/dev/null" % (infasta, musclePath)).read()
            out = getConservedCytosins(
                muscleOut,
                eventObject1['coordinate'],
                eventObject2['coordinate'],
                line[0],
                line[1],
                eventObject1['stand'],
                eventObject2['stand']
            )
            writeFile(sys.argv[4], out)

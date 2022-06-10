'''
Descripttion:
version:
Author: zpliu
Date: 2021-01-11 21:49:13
LastEditors: zpliu
LastEditTime: 2021-01-17 22:19:35
@param:
'''


import os
import re
import sys


def getCIGARItemscode(CIGAR):
    '''
    @Descripttion: get length and CIGAR code
    @param:
    @return:
    '''
    codes = re.split(r'[0-9]+', CIGAR)[1:]
    Length = [int(i) for i in re.split(r'[^0-9]', CIGAR)[0:-1]]
    for code, length in zip(codes, Length):
        yield (code, length)


def getReadCount(BamFile, chromsome, Eventcoortinate):
    '''
    @Descripttion:
    @param:
    @return:
    param {*} BamFile
    param {*} chromsome @str
    param {*} Eventcoortinate @list
    '''
    Coordinate = chromsome+":" + \
        str(Eventcoortinate[2])+"-"+str(Eventcoortinate[3])
    CIGARItems = os.popen(
        'samtools view   -O SAM %s %s |cut -f1,3,4,6' % (BamFile, Coordinate)
    ).read()
    supportIncludeReadCount = 0
    unsupportIncludeReadCount = 0
    otherReads = 0
    # all read
    if not CIGARItems:
        # don't have read
        return (supportIncludeReadCount, unsupportIncludeReadCount, otherReads)
    for item in CIGARItems.split("\n")[0:-1]:
        item = item.split("\t")
        MatchCoordinate = int(item[2])-1
        supportstartSite = 0
        supportEndSite = 0
        unsupportCount = 0
        for code, length in getCIGARItemscode(item[3]):
            if code == 'N' and MatchCoordinate+1+length == Eventcoortinate[0]:
                # match the start sit
                supportstartSite += 1
                MatchCoordinate = MatchCoordinate+length
            elif code == 'N' and MatchCoordinate+1+length == Eventcoortinate[1]:
                # spaning the junction
                supportEndSite += 1
                MatchCoordinate = MatchCoordinate+length
            elif code == 'N' and length == Eventcoortinate[3]-Eventcoortinate[2]-1 and MatchCoordinate+length+1 == Eventcoortinate[3]:
                unsupportCount += 1
                break
            else:
                MatchCoordinate = MatchCoordinate+length
        # decide this read
        if unsupportCount != 0:
            unsupportIncludeReadCount += 1
        elif max((supportstartSite, supportEndSite)) != 0:
            supportIncludeReadCount += 1
        else:
            otherReads += 1
    return (supportIncludeReadCount, unsupportIncludeReadCount, otherReads)


if __name__ == "__main__":
    '''
    @Descripttion:
    @param: inputItem@ str 'Gorai.005G209600;RI:Chr05:59068284-59069242:-'
    @param: RI:chr05:start-1-end:stand
    @return: event  eventLength supportIncludeReadCount unsupportIncludeReadCountotherReads
    '''
    out = []
    EventFile = sys.argv[1]
    BamFile = sys.argv[2]
    outFile = sys.argv[3]
    with open(EventFile, 'r') as File:
        for line in File:
            line = line.strip("\n")
            chromsome = line.split(":")[1]
            supportStartSite = int(line.split(":")[2].split("-")[1])
            supportEndSite = int(line.split(":")[3].split("-")[0])
            unsupportStartSite = int(line.split(":")[2].split("-")[0])
            unsupportEndSite = int(line.split(":")[3].split("-")[1])
            readArray = getReadCount(
                BamFile, chromsome, (supportStartSite, supportEndSite, unsupportStartSite, unsupportEndSite))
            try:
                PSI = readArray[0] / \
                    (readArray[0]+readArray[1])
            except ZeroDivisionError:
                PSI = None
            eventLength = supportEndSite-supportStartSite
            out.append(line+"\t"+str(eventLength) + "\t" +
                       str(PSI)+"\t"+str(readArray[0]+readArray[1])+"\n")
    with open(outFile, 'w') as File:
        for item in out:
            File.write(item)

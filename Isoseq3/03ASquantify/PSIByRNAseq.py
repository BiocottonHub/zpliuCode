'''
Descripttion: 
version: 
Author: zpliu
Date: 2021-01-11 21:49:13
LastEditors: zpliu
LastEditTime: 2021-02-17 13:55:54
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
    def readInEvent(MatchCoordinate, length, Eventcoortinate):

        if MatchCoordinate < Eventcoortinate[0] and MatchCoordinate+length-Eventcoortinate[0] >= 5:
            return True
        elif MatchCoordinate+length > Eventcoortinate[1] and MatchCoordinate - Eventcoortinate[1] <= -5:
            return True
        elif MatchCoordinate > Eventcoortinate[0] and MatchCoordinate+length < Eventcoortinate[1]:
            return True
        else:
            return False

    Coordinate = chromsome+":" + \
        str(Eventcoortinate[0])+"-"+str(Eventcoortinate[1])
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
        startSite = 0
        endSite = 0
        supportCount = 0
        for code, length in getCIGARItemscode(item[3]):
            if code == 'M' and MatchCoordinate+length == Eventcoortinate[0]:
                # match the start sit
                startSite += 1
                MatchCoordinate = MatchCoordinate+length
            elif code == 'N' and MatchCoordinate+length == Eventcoortinate[1]:
                # spaning the junction
                endSite += 1
                MatchCoordinate = MatchCoordinate+length
            elif code == 'M' and readInEvent(MatchCoordinate, length, Eventcoortinate):
                supportCount += 1
                MatchCoordinate = MatchCoordinate+length
            elif code == 'N' and MatchCoordinate < Eventcoortinate[0] and MatchCoordinate+length > Eventcoortinate[1]:
                # the read spaning the intron
                startSite += 1
                endSite += 1
            else:
                MatchCoordinate = MatchCoordinate+length

        # decide this read
        if supportCount != 0:
            supportIncludeReadCount += 1
        elif min((startSite, endSite)) != 0:
            unsupportIncludeReadCount += 1
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
            startSite = int(line.split(":")[2].split(
                "-")[0])-1  # exon and intron junction site
            endSite = int(line.split(":")[2].split("-")[1])
            readArray = getReadCount(
                BamFile, chromsome, (startSite, endSite))
            eventLength = endSite-startSite
            # only apply to RI
            normalsupportReadCount = readArray[0]*150*2/(150*2+eventLength)
            try:
                PSI = normalsupportReadCount / \
                    (normalsupportReadCount+readArray[1])
            except ZeroDivisionError:
                PSI = None
            out.append(line+"\t"+str(eventLength) + "\t" +
                       str(PSI)+"\t"+str(normalsupportReadCount)+"\n")
    with open(outFile, 'w') as File:
        for item in out:
            File.write(item)

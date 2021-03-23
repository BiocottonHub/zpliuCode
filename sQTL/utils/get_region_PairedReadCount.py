'''
Descripttion:
version:
Author: zpliu
Date: 2021-03-19 19:50:51
LastEditors: zpliu
LastEditTime: 2021-03-23 20:59:09
@param:
'''
from utils.parse_bam_flag import parser_flag
import re
import os
import sys


def getCIGARItemscode(CIGAR):
    '''
    @Descripttion: get length and CIGAR code
    @param:
    @return:
    #! '7S8M72844N128M'
    # ? [(7,S),(8,M)...]
    '''

    codes = re.split(r'[0-9]+', CIGAR)[1:]
    Length = [int(i) for i in re.split(r'[^0-9]', CIGAR)[0:-1]]
    for code, length in zip(codes, Length):
        yield (code, length)


def filterRead(FlagStr, MappingQuality):
    flagList = parser_flag(FlagStr)
    # ? filter reads
    '''
    1. mapping quality >20
    2. paired mapping
    3.
    '''
    if int(MappingQuality) < 20:
        return None
    elif 1 not in flagList:
        # not read paired
        return None
    elif 4 in flagList:
        # * unmapping to reference
        return None
    elif 8 in flagList:
        # * mate read unmapping to reference
        return None
    elif 256 in flagList or 1024 in flagList or 2048 in flagList:
        # * other reson to filter the reads
        return None
    elif 16 in flagList and 64 in flagList:
        # * first stand on the reverse stand '-'
        # * mate stand on the  '+'
        return '+'
    elif 32 in flagList and 128 in flagList:
        # * meta stand on the reverse stand '+'
        # * first stand on the '-'
        return '+'
    elif 32 in flagList and 64 in flagList:
        # * meta stand on the reverse stand '-'
        # * first stand on the  '+'
        return '-'
    elif 16 in flagList and 128 in flagList:
        # * meta stand on the  stand '-'
        # * first stand on reverse the  '+'
        return '-'
    else:
        return None


def getReadcount(chromsome, start, end, stand, BamFile, featureId: str):
    '''
    #! get read count
    '''
    # ? samtools view location
    location = chromsome+":"+str(start)+"-"+str(end)
    CIGARItems = os.popen(
        'samtools view   -O BAM %s %s |samtools sort -n -O SAM|cut -f1,2,3,4,5,6,8' % (
            BamFile, location)
    ).read()
    # ? filter the @SQ header
    CIGARItems = [i for i in CIGARItems.split(
        "\n")[0:-1] if re.match('^[^@]', i)]

    if not CIGARItems:
        return (chromsome, start, end, stand, featureId, 0)
    else:
        record = {}
        for item in CIGARItems:
            # ? all reads
            # print(item)
            item = item.split("\t")
            if filterRead(item[1], item[4]) == stand:
                # ? read mapping to the same stand
                #! Continuous same stand but unpaired
                record[item[0]] = record.get(item[0], 0)+1
        # ? paired reads
        readcount = len(list(filter(lambda x: x == 2, record.values())))
        return (chromsome, start, end, stand, featureId, readcount)


if __name__ == "__main__":
    BameFile = sys.argv[1]
    featureId = sys.argv[2]
    # ? location style:
    # ? chr01:111-2222:+
    location = sys.argv[3]
    chromsome = location.split(":")[0]
    start, end = [int(i) for i in location.split(":")[1].split("-")]
    stand = location.split(":")[-1]
    print(getReadcount(
        chromsome,
        start,
        end,
        stand,
        BameFile,
        featureId
    ))

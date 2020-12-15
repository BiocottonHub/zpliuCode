from getorfbylength import getisoformORF
from readgtf import getTranscriptInfo
import argparse
import re


class genes:
    def __init__(self, geneName):
        self.transcriptsORF = []

    def settranscriptORF(self, ORFMessage, ORFAbsoluteCoordinate):
        ORFMessage['AbsoluteStart'] = ORFAbsoluteCoordinate[0]
        ORFMessage['AbsoluteEnd'] = ORFAbsoluteCoordinate[1]
        self.transcriptsORF.append(ORFMessage)

    def compareORF(self, ORFMessage, ORFAbsoluteCoordinate):
        '''
        将PacBio测序的序列与参考基因组序列进行比较
        '''
        # print(self.transcriptsORF)
        # print(ORFAbsoluteCoordinate)
        # print(ORFMessage)
        frameshiftCount = 0
        isoformCount = len(self.transcriptsORF)
        for item in self.transcriptsORF:
            # 编码框发生偏移
            if ORFAbsoluteCoordinate[0] != item['AbsoluteStart'] and ORFAbsoluteCoordinate[1] != item['AbsoluteEnd']:
                frameshiftCount += 1
        if frameshiftCount == isoformCount:
            return {
                'name': ORFMessage['name'],
                'length': ORFMessage['Relativeend']-ORFMessage['Relativestart']+1,
                'start': ORFAbsoluteCoordinate[0],
                'end': ORFAbsoluteCoordinate[1],
                # 粗略的拿最后一个转录本的信息进行比较
                # 'refstart': self.transcriptsORF[0]['AbsoluteStart'],
                # 'refEnd': self.transcriptsORF[0]['AbsoluteEnd'],
                'noreference': 'noreference'
            }
        for item in self.transcriptsORF:
            # 起始位点相同情况
            if item['AbsoluteStart'] == ORFAbsoluteCoordinate[0] and item['AbsoluteEnd'] == ORFAbsoluteCoordinate[1] and ORFMessage['sequence'] == item['sequence']:
                return {
                    'name': ORFMessage['name'],
                    'length': ORFMessage['Relativeend']-ORFMessage['Relativestart']+1,
                    'start': ORFAbsoluteCoordinate[0],
                    'end': ORFAbsoluteCoordinate[1],
                    'ref': item['name'],
                    'refstart': item['AbsoluteStart'],
                    'refEnd': item['AbsoluteEnd'],
                    'frameshift': 'noframeshift',
                    'earlyStop': 'noearlyStop',
                    'laterStop': 'nolaterStop',
                    'earlyStart': 'noearlyStart',
                    'laterStart': 'nolaterStart'

                }
            if item['AbsoluteStart'] == ORFAbsoluteCoordinate[0] and item['AbsoluteEnd'] == ORFAbsoluteCoordinate[1] and ORFMessage['sequence'] != item['sequence']:
              # 起始和终止是一致的，没有造成移码突变，就可能是AS那一段序列也能被翻译出来
                return {
                    'name': ORFMessage['name'],
                    'length': ORFMessage['Relativeend']-ORFMessage['Relativestart']+1,
                    'start': ORFAbsoluteCoordinate[0],
                    'end': ORFAbsoluteCoordinate[1],
                    'ref': item['name'],
                    'refstart': item['AbsoluteStart'],
                    'refEnd': item['AbsoluteEnd'],
                    'inframeChange': 'inframeChange'

                }
            if item['AbsoluteEnd'] > item['AbsoluteStart'] and item['AbsoluteStart'] == ORFAbsoluteCoordinate[0] and item['AbsoluteEnd'] > ORFAbsoluteCoordinate[1]:  # 正链
                if referenceTranscriptInfo[item['name']].getcDNAlengthBycoordinate(ORFAbsoluteCoordinate[0], ORFAbsoluteCoordinate[1], item['AbsoluteStart'], item['AbsoluteEnd'], PacBioTranscriptInfo[ORFMessage['name']].getExonCoordinate()) % 3 == 0:
                    return {
                        'name': ORFMessage['name'],
                        'length': ORFMessage['Relativeend']-ORFMessage['Relativestart']+1,
                        'start': ORFAbsoluteCoordinate[0],
                        'end': ORFAbsoluteCoordinate[1],
                        'ref': item['name'],
                        'refstart': item['AbsoluteStart'],
                        'refEnd': item['AbsoluteEnd'],
                        'frameshift': 'noframeshift',
                        'earlyStop': 'earlyStop',
                        'laterStop': 'nolaterStop',
                        'earlyStart': 'noearlyStart',
                        'laterStart': 'nolaterStart'
                    }
                else:
                    return {
                        'name': ORFMessage['name'],
                        'length': ORFMessage['Relativeend']-ORFMessage['Relativestart']+1,
                        'start': ORFAbsoluteCoordinate[0],
                        'end': ORFAbsoluteCoordinate[1],
                        'ref': item['name'],
                        'refstart': item['AbsoluteStart'],
                        'refEnd': item['AbsoluteEnd'],
                        'frameshift': 'frameshift',
                        'earlyStop': 'earlyStop',
                        'laterStop': 'nolaterStop',
                        'earlyStart': 'noearlyStart',
                        'laterStart': 'nolaterStart'
                    }
            if item['AbsoluteEnd'] < item['AbsoluteStart'] and item['AbsoluteStart'] == ORFAbsoluteCoordinate[0] and item['AbsoluteEnd'] > ORFAbsoluteCoordinate[1]:  # 负链
                if referenceTranscriptInfo[item['name']].getcDNAlengthBycoordinate(ORFAbsoluteCoordinate[0], ORFAbsoluteCoordinate[1], item['AbsoluteStart'], item['AbsoluteEnd'], PacBioTranscriptInfo[ORFMessage['name']].getExonCoordinate()) % 3 == 0:
                    return {
                        'name': ORFMessage['name'],
                        'length': ORFMessage['Relativeend']-ORFMessage['Relativestart']+1,
                        'start': ORFAbsoluteCoordinate[0],
                        'end': ORFAbsoluteCoordinate[1],
                        'ref': item['name'],
                        'refstart': item['AbsoluteStart'],
                        'refEnd': item['AbsoluteEnd'],
                        'frameshift': 'noframeshift',
                        'earlyStop': 'noearlyStop',
                        'laterStop': 'laterStop',
                        'earlyStart': 'noearlyStart',
                        'laterStart': 'nolaterStart'
                    }
                else:
                    return {
                        'name': ORFMessage['name'],
                        'length': ORFMessage['Relativeend']-ORFMessage['Relativestart']+1,
                        'start': ORFAbsoluteCoordinate[0],
                        'end': ORFAbsoluteCoordinate[1],
                        'ref': item['name'],
                        'refstart': item['AbsoluteStart'],
                        'refEnd': item['AbsoluteEnd'],
                        'frameshift': 'frameshift',
                        'earlyStop': 'noearlyStop',
                        'laterStop': 'laterStop',
                        'earlyStart': 'noearlyStart',
                        'laterStart': 'nolaterStart'
                    }
            if item['AbsoluteEnd'] > item['AbsoluteStart'] and item['AbsoluteStart'] == ORFAbsoluteCoordinate[0] and item['AbsoluteEnd'] < ORFAbsoluteCoordinate[1]:
                if referenceTranscriptInfo[item['name']].getcDNAlengthBycoordinate(ORFAbsoluteCoordinate[0], ORFAbsoluteCoordinate[1], item['AbsoluteStart'], item['AbsoluteEnd'], PacBioTranscriptInfo[ORFMessage['name']].getExonCoordinate()) % 3 == 0:
                    return {
                        'name': ORFMessage['name'],
                        'length': ORFMessage['Relativeend']-ORFMessage['Relativestart']+1,
                        'start': ORFAbsoluteCoordinate[0],
                        'end': ORFAbsoluteCoordinate[1],
                        'ref': item['name'],
                        'refstart': item['AbsoluteStart'],
                        'refEnd': item['AbsoluteEnd'],
                        'frameshift': 'noframeshift',
                        'earlyStop': 'noearlyStop',
                        'laterStop': 'laterStop',
                        'earlyStart': 'noearltStart',
                        'laterStart': 'nolaterStart'
                    }
                else:
                    return {
                        'name': ORFMessage['name'],
                        'length': ORFMessage['Relativeend']-ORFMessage['Relativestart']+1,
                        'start': ORFAbsoluteCoordinate[0],
                        'end': ORFAbsoluteCoordinate[1],
                        'ref': item['name'],
                        'refstart': item['AbsoluteStart'],
                        'refEnd': item['AbsoluteEnd'],
                        'frameshift': 'frameshift',
                        'earlyStop': 'noearlyStop',
                        'laterStop': 'laterStop',
                        'earlyStart': 'noearltStart',
                        'laterStart': 'nolaterStart'
                    }
            if item['AbsoluteEnd'] < item['AbsoluteStart'] and item['AbsoluteStart'] == ORFAbsoluteCoordinate[0] and item['AbsoluteEnd'] < ORFAbsoluteCoordinate[1]:
                if referenceTranscriptInfo[item['name']].getcDNAlengthBycoordinate(ORFAbsoluteCoordinate[0], ORFAbsoluteCoordinate[1], item['AbsoluteStart'], item['AbsoluteEnd'], PacBioTranscriptInfo[ORFMessage['name']].getExonCoordinate()) % 3 == 0:
                    return {
                        'name': ORFMessage['name'],
                        'length': ORFMessage['Relativeend']-ORFMessage['Relativestart']+1,
                        'start': ORFAbsoluteCoordinate[0],
                        'end': ORFAbsoluteCoordinate[1],
                        'ref': item['name'],
                        'refstart': item['AbsoluteStart'],
                        'refEnd': item['AbsoluteEnd'],
                        'frameshift': 'noframeshift',
                        'earlyStop': 'earlyStop',
                        'laterStop': 'nolaterStop',
                        'earlyStart': 'noearlyStart',
                        'laterStart': 'nolaterStart'
                    }
                else:
                    return {
                        'name': ORFMessage['name'],
                        'length': ORFMessage['Relativeend']-ORFMessage['Relativestart']+1,
                        'start': ORFAbsoluteCoordinate[0],
                        'end': ORFAbsoluteCoordinate[1],
                        'ref': item['name'],
                        'refstart': item['AbsoluteStart'],
                        'refEnd': item['AbsoluteEnd'],
                        'frameshift': 'frameshift',
                        'earlyStop': 'earlyStop',
                        'laterStop': 'nolaterStop',
                        'earlyStart': 'noearlyStart',
                        'laterStart': 'nolaterStart'
                    }
            # 终止位点相同情况
            if item['AbsoluteEnd'] > item['AbsoluteStart'] and item['AbsoluteEnd'] == ORFAbsoluteCoordinate[1] and item['AbsoluteStart'] > ORFAbsoluteCoordinate[0]:
                if referenceTranscriptInfo[item['name']].getcDNAlengthBycoordinate(ORFAbsoluteCoordinate[0], ORFAbsoluteCoordinate[1], item['AbsoluteStart'], item['AbsoluteEnd'], PacBioTranscriptInfo[ORFMessage['name']].getExonCoordinate()) % 3 == 0:
                    return {
                        'name': ORFMessage['name'],
                        'length': ORFMessage['Relativeend']-ORFMessage['Relativestart']+1,
                        'start': ORFAbsoluteCoordinate[0],
                        'end': ORFAbsoluteCoordinate[1],
                        'ref': item['name'],
                        'refstart': item['AbsoluteStart'],
                        'refEnd': item['AbsoluteEnd'],
                        'frameshift': 'noframeshift',
                        'earlyStop': 'noearlyStop',
                        'laterStop': 'nolaterStop',
                        'earlyStart': 'earlyStart',
                        'laterStart': 'nolaterStart'
                    }
                else:
                    return {
                        'name': ORFMessage['name'],
                        'length': ORFMessage['Relativeend']-ORFMessage['Relativestart']+1,
                        'start': ORFAbsoluteCoordinate[0],
                        'end': ORFAbsoluteCoordinate[1],
                        'ref': item['name'],
                        'refstart': item['AbsoluteStart'],
                        'refEnd': item['AbsoluteEnd'],
                        'frameshift': 'frameshift',
                        'earlyStop': 'noearlyStop',
                        'laterStop': 'nolaterStop',
                        'earlyStart': 'earlyStart',
                        'laterStart': 'nolaterStart'
                    }
            if item['AbsoluteEnd'] > item['AbsoluteStart'] and item['AbsoluteEnd'] == ORFAbsoluteCoordinate[1] and item['AbsoluteStart'] < ORFAbsoluteCoordinate[0]:
                if referenceTranscriptInfo[item['name']].getcDNAlengthBycoordinate(ORFAbsoluteCoordinate[0], ORFAbsoluteCoordinate[1], item['AbsoluteStart'], item['AbsoluteEnd'], PacBioTranscriptInfo[ORFMessage['name']].getExonCoordinate()) % 3 == 0:
                    return {
                        'name': ORFMessage['name'],
                        'length': ORFMessage['Relativeend']-ORFMessage['Relativestart']+1,
                        'start': ORFAbsoluteCoordinate[0],
                        'end': ORFAbsoluteCoordinate[1],
                        'ref': item['name'],
                        'refstart': item['AbsoluteStart'],
                        'refEnd': item['AbsoluteEnd'],
                        'frameshift': 'noframeshift',
                        'earlyStop': 'noearlyStop',
                        'laterStop': 'nolaterStop',
                        'earlyStart': 'noearlyStart',
                        'laterStart': 'laterStart'
                    }
                else:
                    return {
                        'name': ORFMessage['name'],
                        'length': ORFMessage['Relativeend']-ORFMessage['Relativestart']+1,
                        'start': ORFAbsoluteCoordinate[0],
                        'end': ORFAbsoluteCoordinate[1],
                        'ref': item['name'],
                        'refstart': item['AbsoluteStart'],
                        'refEnd': item['AbsoluteEnd'],
                        'frameshift': 'frameshift',
                        'earlyStop': 'noearlyStop',
                        'laterStop': 'nolaterStop',
                        'earlyStart': 'noearlyStart',
                        'laterStart': 'laterStart'
                    }
            if item['AbsoluteEnd'] < item['AbsoluteStart'] and item['AbsoluteEnd'] == ORFAbsoluteCoordinate[1] and item['AbsoluteStart'] > ORFAbsoluteCoordinate[0]:
                if referenceTranscriptInfo[item['name']].getcDNAlengthBycoordinate(ORFAbsoluteCoordinate[0], ORFAbsoluteCoordinate[1], item['AbsoluteStart'], item['AbsoluteEnd'], PacBioTranscriptInfo[ORFMessage['name']].getExonCoordinate()) % 3 == 0:
                    return {
                        'name': ORFMessage['name'],
                        'length': ORFMessage['Relativeend']-ORFMessage['Relativestart']+1,
                        'start': ORFAbsoluteCoordinate[0],
                        'end': ORFAbsoluteCoordinate[1],
                        'ref': item['name'],
                        'refstart': item['AbsoluteStart'],
                        'refEnd': item['AbsoluteEnd'],
                        'frameshift': 'noframeshift',
                        'earlyStop': 'noearlyStop',
                        'laterStop': 'nolaterStop',
                        'earlyStart': 'noearlyStart',
                        'laterStart': 'laterStart'
                    }
                else:
                    return {
                        'name': ORFMessage['name'],
                        'length': ORFMessage['Relativeend']-ORFMessage['Relativestart']+1,
                        'start': ORFAbsoluteCoordinate[0],
                        'end': ORFAbsoluteCoordinate[1],
                        'ref': item['name'],
                        'refstart': item['AbsoluteStart'],
                        'refEnd': item['AbsoluteEnd'],
                        'frameshift': 'frameshift',
                        'earlyStop': 'noearlyStop',
                        'laterStop': 'nolaterStop',
                        'earlyStart': 'noearlyStart',
                        'laterStart': 'laterStart'
                    }

            if item['AbsoluteEnd'] > item['AbsoluteStart'] and item['AbsoluteEnd'] == ORFAbsoluteCoordinate[1] and item['AbsoluteStart'] < ORFAbsoluteCoordinate[0]:
                if referenceTranscriptInfo[item['name']].getcDNAlengthBycoordinate(ORFAbsoluteCoordinate[0], ORFAbsoluteCoordinate[1], item['AbsoluteStart'], item['AbsoluteEnd'], PacBioTranscriptInfo[ORFMessage['name']].getExonCoordinate()) % 3 == 0:
                    return {
                        'name': ORFMessage['name'],
                        'length': ORFMessage['Relativeend']-ORFMessage['Relativestart']+1,
                        'start': ORFAbsoluteCoordinate[0],
                        'end': ORFAbsoluteCoordinate[1],
                        'ref': item['name'],
                        'refstart': item['AbsoluteStart'],
                        'refEnd': item['AbsoluteEnd'],
                        'frameshift': 'noframeshift',
                        'earlyStop': 'noearlyStop',
                        'laterStop': 'nolaterStop',
                        'earlyStart': 'earlyStart',
                        'laterStart': 'nolaterStart'
                    }
                else:
                    return {
                        'name': ORFMessage['name'],
                        'length': ORFMessage['Relativeend']-ORFMessage['Relativestart']+1,
                        'start': ORFAbsoluteCoordinate[0],
                        'end': ORFAbsoluteCoordinate[1],
                        'ref': item['name'],
                        'refstart': item['AbsoluteStart'],
                        'refEnd': item['AbsoluteEnd'],
                        'frameshift': 'frameshift',
                        'earlyStop': 'noearlyStop',
                        'laterStop': 'nolaterStop',
                        'earlyStart': 'earlyStart',
                        'laterStart': 'nolaterStart'
                    }
        return {'name': ORFMessage['name'], 'length': ORFMessage['Relativeend']-ORFMessage['Relativestart']+1, }


'''
首先根据参考基因组的注释信息，
将同一个基因预测的ORF信息存成一个对象；
遍历每个PacBio预测的ORF信息，将它与参考基因的ORF信息进行比较存在以下几种情况
+ 氨基酸序列与参考基因组中的某个转录本一致
+ 起始的位置相同，但是相比于参考基因组转录本提前终止，这就是frameshift发生偏移了
+ 起始位置相同，但是相比于参考基因组转录本
##判断是否发生frameshift
如果起始位点相同，根据PacBio终止位点的绝对坐标，提取参考基因组对应区域的cDNA长度;将这个长度与PacBio的ORF长度进行比较；是否是相差3的倍数；
终止位点相同的情况类似进行分析
'''
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="")
    parser.add_argument("-refgtf", help="reference gtf file")
    parser.add_argument("-PBgtf", help="PacBio gtf file")
    parser.add_argument("-reforf", help="reference ORF")
    parser.add_argument("-PBorf", help="PacBio ORF")
    parser.add_argument("-o", help="out put File")
    args = parser.parse_args()
    # 转录本的ORF信息
    PacBioIsofor = getisoformORF(args.PBgtf, args.PBorf)
    referIsofor = getisoformORF(args.refgtf, args.reforf)
    referenceTranscriptInfo = getTranscriptInfo(args.refgtf)
    PacBioTranscriptInfo = getTranscriptInfo(args.PBgtf)
    genePattern = r'gene_id "([^"]*)"'
    transcriptPattern = r'transcript_id "([^"]*)"'
    tmp = ''
    geneDict = {}
    tmp2 = []
    with open(args.refgtf, 'r') as File:
        for line in File.readlines():
            line = line.strip("\n").split("\t")
            geneId = re.search(genePattern, line[8]).group(1)
            transcriptId = re.search(transcriptPattern, line[8]).group(1)
            if transcriptId == tmp:  # 跳过多个外显子区域
                continue
            else:
                tmp = transcriptId
            if geneId not in tmp2:
                geneDict[geneId] = genes(geneId)
                tmp2.append(geneId)  # 标记这个基因已经初始化了
            try:
                geneDict[geneId].settranscriptORF(
                    referIsofor[transcriptId].getORF(), referIsofor[transcriptId].getAbsoulteORF())
            except KeyError:
                # 基因的转录本没有预测出ORF
                pass
    # 开始与参考基因组的注释进行比较
    # print(geneDict['Ghir_A01G009300'].compareORF(
    #     PacBioIsofor['PB.5777.1'].getORF(), PacBioIsofor['PB.5777.1'].getAbsoulteORF()))
    with open(args.o, 'w') as File:
        for key in PacBioIsofor.keys():
            geneName = PacBioIsofor[key].getgeneName()
            OutDict = geneDict[geneName].compareORF(
                PacBioIsofor[key].getORF(), PacBioIsofor[key].getAbsoulteORF())
            # print(OutDict)
            File.write(geneName+"\t"+"\t".join([str(OutDict[i])
                                                for i in OutDict.keys()])+"\n")

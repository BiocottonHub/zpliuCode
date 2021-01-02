'''
Descripttion:
version:
Author: zpliu
Date: 2020-12-29 09:45:06
LastEditors: zpliu
LastEditTime: 2021-01-02 20:55:31
@param:
'''
import re
import os
import argparse
from multiprocessing import Process, Queue
'''
function:
@getEventLength: get AS event length
@readHomolog: get homolog gene dict
@readFastaFile: get gene fasta sequence
@getKmerSequence: get kmer sequence from gene sequence
@readgeneCoordinate: get gene coordinate 
@getASsequence: get AS flagment sequence
@muscle: run muscle form two sequence
@runASevent: get out form a AS list
'''


class multiProcessASKmer(Process):
    def __init__(self, processId: int, queue: Queue, ASevents: list, homologDict: dict, querygeneCoordinateDict: dict, datageneCoordinateDict: dict, homologGeneFasta: dict):
        # inherit form process object
        super(multiProcessASKmer, self).__init__()
        self.__id = processId
        self.__ASevents = ASevents
        self.__homologDict = homologDict
        self.__querygeneCoordinateDict = querygeneCoordinateDict
        self.__datageneCoordinateDict = datageneCoordinateDict
        self.__homologGeneFasta = homologGeneFasta
        self.__queue = queue  # single process out data

    def run(self):
        print("Process: " + str(self.__id) +
              ": begin to deal with ASs...")
        # get each process out
        outlist = runASevent(self.__ASevents, self.__homologDict, self.__querygeneCoordinateDict,
                             self.__datageneCoordinateDict, self.__homologGeneFasta)
        for item in outlist:
            # put data into Queue
            self.__queue.put(item, block=True, timeout=None)
        print("Process: " + str(self.__id) +
              ": finished...")


def getEventLength(ASevent):
    '''
    @Descripttion:  according AS file to get AS coordinate
    @param: ASevent @str
    @return: AS message @dict
    '''
    if re.search('RI', ASevent):
        local = re.search(r':([0-9]*)-([0-9]*):', ASevent)

        return {
            'length': abs(int(local.group(1))-int(local.group(2))),
            'coordinate': [int(local.group(1)), int(local.group(2))]
        }
    if re.search('SE', ASevent):
        local = re.search(r'-([0-9]*):([0-9]*)-', ASevent)
        return {
            'length': abs(int(local.group(1))-int(local.group(2))),
            'coordinate': [int(local.group(1)), int(local.group(2))]
        }
    if re.search('A3', ASevent):
        local = re.search(r':([0-9]*)-([0-9]*):([0-9]*)-([0-9]*):', ASevent)
        tmp = sorted([int(local.group(1)), int(local.group(2)),
                      int(local.group(3)), int(local.group(4))])
        return {
            'length': tmp[-1]-tmp[0],
            'coordinate': [tmp[0], tmp[-1]]
        }
    if re.search('A5', ASevent):
        local = re.search(r':([0-9]*)-([0-9]*):([0-9]*)-([0-9]*):', ASevent)
        tmp = sorted([int(local.group(1)), int(local.group(2)),
                      int(local.group(3)), int(local.group(4))])
        return {
            'length': tmp[-1]-tmp[0],
            'coordinate': [tmp[0], tmp[-1]]
        }


def readHomolog(homologGeneFile):
    '''
    @Descripttion: get homologous gene
    @param: homologGeneFile @dict
    @return: homolog @dict
    '''
    homologGene = {}
    with open(homologGeneFile, 'r') as File:
        for line in File:
            line = line.strip("\n").split("\t")
            homologGene[line[2]] = line[3]
            homologGene[line[3]] = line[2]
    return homologGene


def readFastaFile(geneFile):
    '''
    @Descripttion: get gene sequence from fasta file
    @param: gene fasta file @str
    @return: gene sequence @dict
    '''
    out = {}
    with open(geneFile, 'r') as File:
        lines = File.readlines()
        for i in range(0, len(lines), 2):
            out[lines[i].strip("\n").strip(">")] = lines[i+1].strip("\n")
    return out


def getKmerSequence(sequence, geneId, k, geneCoordinate):
    '''
    @Descripttion: from gene sequence to get kmer sequence
    @param: gene fasta sequence @str
    @param: geneId @str
    @param: AS length as k-mer length  @int
    @param: geneCoordinate Message contained  @dict
    @return:  kmer sequence @dict
    '''
    end = len(sequence)
    out = {}
    for i in range(0, (end+1)-k+1):
        seq = sequence[i:i+k]
        if geneCoordinate['stand'] == "+":
            name = ">"+geneId+";" + \
                geneCoordinate['chr'] + ":" + \
                str(geneCoordinate['start']+i)+"-" + \
                str(geneCoordinate['start']+i+k-1)
        else:
            # - with the absolute kmer
            name = ">"+geneId+";" + \
                geneCoordinate['chr'] + ":" + \
                str(geneCoordinate['end']-i-k)+"-" + \
                str(geneCoordinate['end']-i)
        out[name] = name+"\n"+seq
    return out


def readgeneCoordinate(geneCoordiante):
    '''
    @Descripttion: get gene Message
    @param: geneCoordiante @str
    @return: @dict{
      'geneID'{
        'chr' @str
        'start' @int
        'end' @int
      }
    }
    '''
    out = {}
    with open(geneCoordiante, 'r') as File:
        for line in File:
            line = line.strip("\n").split("\t")
            out[line[3]] = {
                'chr': line[0],
                'start': int(line[1]),
                'end': int(line[2]),
                'stand': line[-1]
            }
    return out


def getASsequence(ASeventName, AScoordinate, genesequence, geneCoordinate):
    '''
    @Descripttion: According gene sequence to get AS sequence
    @param:
    @return:
    '''
    start1, end1 = sorted(AScoordinate)
    # AS enevt out range of gene
    if start1 < geneCoordinate['start'] or end1 > geneCoordinate['end']:
        return None
    else:
        if geneCoordinate['stand'] == "+":
            start = start1-geneCoordinate['start']+1
            end = start+end1-start1
            return ">"+ASeventName+"\n"+genesequence[start:end]
        else:
            start = geneCoordinate['end']-end1
            end = start+end1-start1
            return ">"+ASeventName+"\n"+genesequence[start:end]


def muscle(ASsequence, kmerSequence):
    infasta = ASsequence+"\n"+kmerSequence+"\n"
    out = os.popen(
        "printf \"%s\"|~/software/muscle3.8.31_i86linux64  -clw  2>/dev/null" % infasta).read()
    return len(re.findall(r'\*', out))


def runASevent(ASevents, homologDict, querygeneCoordinateDict, datageneCoordinateDict, homologGeneFasta):
    '''
    @Descripttion: get the most identity AS kmer by AS list
    @param: ASevents @list
    @param: homolog gene dict @dict
    @param: querygene AS gene coordiante @dict
    @param: homolog gene coordinate @dict
    @param: homolog gene fasta sequence @dict
    @return: AS coordinate and kmer coordinate and identity nucl @list
    '''
    out = []
    for line in ASevents:
        line = line.strip("\n").split("\t")[2]
        ASevent = getEventLength(line)
        ASgeneId = re.split(";", line)[0]
        searchGeneId = homologDict[ASgeneId]
        ASLength = ASevent['length']
        ASCoordinate = ASevent['coordinate']
        geneCoordinate = querygeneCoordinateDict[ASgeneId]
        ASsequence = getASsequence(
            line, ASCoordinate, homologGeneFasta[ASgeneId], geneCoordinate)
        # print(ASsequence)
        if not ASsequence:
            # AS enevt out range of gene
            print("out range of gene:"+line)
            pass
        else:
            searchKmersequence = getKmerSequence(
                homologGeneFasta[searchGeneId], searchGeneId, ASLength, datageneCoordinateDict[searchGeneId])
            max_rank = 0
            max_kmer = ''
            for key, values in searchKmersequence.items():
                tmp = muscle(ASsequence, values)
                if tmp == ASLength:
                    # the most identity kmer
                    max_kmer = key
                    max_rank = tmp
                    break
                elif tmp > max_rank:
                    # iterator other k-mer
                    max_kmer = key
                    max_rank = tmp
                else:
                    pass
            tmp = line+"\t"+max_kmer+"\t"+str(max_rank)+"\t"+str(ASLength)+"\n"
            out.append(tmp)
    return out


def ASmuscle(queryASFile, geneSequenceFile, querygeneCoordinateFile, datageneCoordinateFile, homologFile, processNum: int):
    '''
    @Descripttion: multiple process to deal with AS events
    @return:
    '''
    homologGeneFasta = readFastaFile(geneSequenceFile)
    homologDict = readHomolog(homologFile)
    querygeneCoordinateDict = readgeneCoordinate(querygeneCoordinateFile)
    datageneCoordinateDict = readgeneCoordinate(datageneCoordinateFile)
    out = []
    processList = []
    workQueue = Queue()
    with open(queryASFile, 'r') as File:
        ASevents = File.readlines()
    average = int(len(ASevents)/processNum)
    for processId in range(0, processNum):
        if processId == processNum-1:
            start = processId*average
            end = len(ASevents)
        else:
            start = processId*average
            end = (processId+1)*average
        myprocess = multiProcessASKmer(
            processId, workQueue, ASevents[start:end], homologDict, querygeneCoordinateDict, datageneCoordinateDict, homologGeneFasta)
        myprocess.start()
        processList.append(myprocess)
    for process in processList:
        process.join()
    while True:
        try:
            # get child process Data
            out.append(workQueue.get(block=True, timeout=1))
        except:
            break
    # out += runASevent(ASevents, homologDict, querygeneCoordinateDict, datageneCoordinateDict, homologGeneFasta)
    return out


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-AS', help="query AS File")
    parser.add_argument('-querygene', help="query gene coordinate")
    parser.add_argument('-datagene', help="data gene coordinate")
    parser.add_argument('-genefasta', help="data gene fasta")
    parser.add_argument('-homolog', help="homolog Gene File")
    parser.add_argument('-p', help="process number")
    parser.add_argument('-out', help="out put File")
    namespace = parser.parse_args()
    out = ASmuscle(namespace.AS, namespace.genefasta,
                   namespace.querygene, namespace.datagene, namespace.homolog, int(namespace.p))
    with open(namespace.out, 'w') as File:
        for item in out:
            File.write(item)

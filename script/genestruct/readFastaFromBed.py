'''
/*
 * @Author: zpliu 
 * @Date: 2020-04-21 14:15:58 
 * @Last Modified by: zpliu
 * @Last Modified time: 2020-04-21 15:37:21
 */
 BED File:
    chrosome start end chain geneId
 fasta file:
    >chr01
    ATATATA
'''
import sys
import re
import pybedtools
#from tqdm import tqdm


def readFastaFile(genomeFile):
    with open(genomeFile, 'r') as File:
        '''
        @Descripttion: get fasta dict
        @param: 
        @return: 
        '''
        genes_dict = {}
        geneId = ''
        for line in File:
            # print(line)
            if re.match(r'^>', line):
                geneId = re.split(r'\s+', line.strip(">").strip("\n"))[0]
                genes_dict[geneId] = genes_dict.get(geneId, '')
            else:
                genes_dict[geneId] += line.strip("\n")
    return genes_dict


def getFastaBypybedtools(Bedstring, FastaFile):
    '''
    @Descripttion: extrract sequence by Bedtools
    @param: 
    @return: 
    '''
    a = pybedtools.BedTool(Bedstring, from_string=True)
    fasta = pybedtools.example_filename(FastaFile)
    a = a.sequence(fi=fasta, name=True, s=True)
    return open(a.seqfn).read()


def reverseSequence(seq):
    dict1 = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C',
             'N': 'N', 'a': 'T', 't': 'A', 'g': 'C', 'c': 'G'}
    rever_seq1 = [dict1[k] for k in seq[::-1]]
    return(''.join(rever_seq1))


if __name__ == "__main__":
    print(readFastaFile(sys.argv[1]).keys()[0:13])

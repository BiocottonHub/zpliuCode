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
from tqdm import tqdm


def readFastaFile(genomeFile):
    with open(genomeFile, 'r') as File:
        genes_fastas = File.readlines()
        genes_dict = {}
    for i in tqdm(range(0, len(genes_fastas)), desc="read genome file"):
        if re.match("^>", genes_fastas[i]):
            genes_fastas[i] = genes_fastas[i].strip("\n").strip(">")
            index1 = i+1
            sequence = ""
            while re.match("^[^>]", genes_fastas[index1]):
                sequence += genes_fastas[index1].strip("\n")
                index1 += 1
                if index1 == len(genes_fastas):
                    break
            genes_dict[genes_fastas[i]] = sequence
        else:
            continue
    return genes_dict


def reverseSequence(seq):
    dict1 = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C',
             'N': 'N', 'a': 'T', 't': 'A', 'g': 'C', 'c': 'G'}
    rever_seq1 = [dict1[k] for k in seq[::-1]]
    return(''.join(rever_seq1))

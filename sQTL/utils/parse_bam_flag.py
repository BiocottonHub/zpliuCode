'''
Descripttion: 
version: 
Author: zpliu
Date: 2021-03-19 19:37:33
LastEditors: zpliu
LastEditTime: 2021-03-19 20:57:39
@param: 
'''
'''
# ? @abstract the read is paired in sequencing, no matter whether it is mapped in a pair * /
# define BAM_FPAIRED        1
#? @abstract the read is mapped in a proper pair * /
# define BAM_FPROPER_PAIR   2
#? @abstract the read itself is unmapped
conflictive with BAM_FPROPER_PAIR * /
# define BAM_FUNMAP         4
#? @abstract the mate is unmapped * /
# define BAM_FMUNMAP        8
#? @abstract the read is mapped to the reverse strand * /
# define BAM_FREVERSE      16
#? @abstract the mate is mapped to the reverse strand * /
# define BAM_FMREVERSE     32
#? @abstract this is read1 * /
# define BAM_FREAD1        64
#? @abstract this is read2 * /
# define BAM_FREAD2       128
#? @abstract not primary alignment * /
# define BAM_FSECONDARY   256
#? @abstract QC failure * /
# define BAM_FQCFAIL      512
#? @abstract optical or PCR duplicate * /
# define BAM_FDUP        1024
#? @abstract supplementary alignment * /
# define BAM_FSUPPLEMENTARY 2048
'''




import sys
def parser_flag(flagStr: str):
    flagindex = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048]
    BinaryCode = bin(int(flagStr))
    out = []
    for index, value in enumerate(str(BinaryCode)[::-1]):
        if value == "1":
            out.append(flagindex[index])
    return out


if __name__ == "__main__":
    print(parser_flag(sys.argv[1]))

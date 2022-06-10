'''
Descripttion: some usefull function
version: 
Author: zpliu
Date: 2020-11-06 10:33:50
LastEditors: zpliu
LastEditTime: 2020-11-06 10:36:09
'''


def getAllFiles(TopPath, ext):
    '''
    @Descripttion: itera the suffix files in the special directory
    @param: TopPath: @str the full path of the directory
    @param: ext: @str  the suffix of files
    @return: @list[str] the lisf of file's path 
    '''
    filePath = []
    for parentDir, dirnames, filenames in os.walk(TopPath):
        for dirname in dirnames:
            filePath += [os.path.join(os.path.join(parentDir, dirname, i)) for i in filter(
                lambda filename: os.path.splitext(filename)[1] == ext,
                os.listdir(os.path.join(parentDir, dirname))
            )]
    return filePath

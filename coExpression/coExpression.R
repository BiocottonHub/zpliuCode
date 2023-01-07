library(WGCNA)
library(dplyr)
library(glue)
options(stringsAsFactors = FALSE)
#* 读取筛选后的表达谱
args <- commandArgs(trailingOnly = TRUE)
exprssRData=args[1]
power=as.numeric(args[2])
outPrefix=args[3]
nthreads=as.numeric(args[4])

print(glue(
"
#############################################
Start run WGCNA...
#############################################
inputFile: {exprssRData}
select power threshold: {power}
Output File:
    + {outPrefix}-TOM.RData
    + {outPrefix}-net.RData
"
))

#? 变量为datExpre_filter
load(file = exprssRData)


#* 开启多线程
#* 进行共表达分析
enableWGCNAThreads() 
net=blockwiseModules(
    datExpr = datExpre_filter,
    power = power,
    TOMType = "unsigned",
    minModuleSize = 30,
    reassignThreshold = 0,
    mergeCutHeight = 0.25,
    numericLabels = TRUE,
    pamRespectsDendro = FALSE,
    saveTOMs = TRUE,
    saveTOMFileBase = paste(outPrefix,"-TOM",sep=""),
    verbose = 3,
    nThreads=nthreads,
    maxBlockSize=12000
)
#* 保存net文件
save(net,file=paste(outPrefix,"-net.RData",sep=""))
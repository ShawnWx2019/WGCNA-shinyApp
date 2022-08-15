#####################################
# 	Desc: WGCNA for server pipeline
# 	Date: Aug 11, 2022
#   File: Step2. network construction
#	  Author: Shawn Wang
#####################################

# depends -----------------------------------------------------------------

suppressMessages(suppressWarnings(library(getopt)));
suppressMessages(suppressWarnings(library(logr)));
load("step1.RData")
# parameters --------------------------------------------------------------
tmp <- file.path("Step2.log")

lf <- log_open(tmp)

command=matrix(c(
  'help', 'h', 0, 'logic', 'help information',
  'power', 'p' , 1, 'double', 'SFT power.',
  'minModuleSize', 'm' , 1, 'double', 'Minimal module size: The gene number of minimal module. Default: 30',
  'mergeCutHeight', 'c' , 1, 'double', 'Merge cuttree height: tree height lower than this value will be merged. Default: 0.25',
  'maxBlocksize', 'b' , 1, 'double', 'max block size: For block-wised network construction method, the block size. Default: all filtered genes.'
),byrow = T, ncol = 5)
args = getopt(command)

log_print(
  "<== Check parameter ==>"
)

## help information
if (!is.null(args$help)) {
  cat(paste(getopt(command, usage = T), "\n"))
  q(status=1)
}

# default
if (is.null(args$power)){
  args$power = step1_sft$sft
}
power = args$power

if (is.null(args$minModuleSize)){
  args$minModuleSize = 30
}
minModuleSize = args$minModuleSize

if (is.null(args$mergeCutHeight)){
  args$mergeCutHeight = 0.25
}
mergeCutHeight = args$mergeCutHeight

if (is.null(args$maxBlocksize)){
  args$maxBlocksize = ncol(datExpr)
}
maxBlocksize = args$maxBlocksize

# run ---------------------------------------------------------------------
suppressMessages(suppressWarnings(library(WGCNA)));
suppressMessages(suppressWarnings(library(ShinyWGCNA)));
suppressMessages(suppressWarnings(library(tidyverse)));
suppressMessages(suppressWarnings(library(ggprism)));
suppressMessages(suppressWarnings(library(patchwork)));
enableWGCNAThreads(nThreads = threads)
# test

# datExpr <- read.table("~/SynologyDrive/Project/08.WheatLifeCycle/03.Progress/06.transcriptomics/AK_CK_expmat_filtered.xls",header = T,sep = "\t")
# ngenes = nrow(datExpr)
# RcCutoff = 0
# samplePerc = 0
# GeneNum = 50000
# GeneNumCut = 1-(GeneNum/ngenes)
# cutmethod = "MAD"
# rscut = 0.8
# datatype = 'count'
# method = "varianceStabilizingTransformation"

log_print(
  "<== Start! ==>"
)

tryCatch({
  step2_network = getnetwork(
    datExpr = datExpr,
    power = power,
    minModuleSize = minModuleSize,
    mergeCutHeight = mergeCutHeight,
    maxBlocksize = maxBlocksize
  )
  log_print(
    "Success!"
  )
},error = function(e) {
  log_print(
    "Failed!"
  )
  q(status = 1)
})


params = data.frame(
  param = c("datatype","method","RcCutoff","samplePerc","GeneNum","cutmethod","rscut","threads","power","minModuleSize","mergeCutHeight","maxBlocksize"),
  value = c(datatype,method,RcCutoff,samplePerc,GeneNum,cutmethod,rscut,threads,power,minModuleSize,mergeCutHeight,maxBlocksize)
)

log_print("Paramters used in this step: \n")
log_print(params)

save.image("step2.RData")

log_print(
  "< == All work finish ==>"
)

log_close()

writeLines(readLines(lf))

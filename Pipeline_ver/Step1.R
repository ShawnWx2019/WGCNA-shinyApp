#####################################
# 	Desc: WGCNA for server pipeline
# 	Date: Aug 11, 2022
#   File: Data cleaning and get sft
#   Author: Shawn Wang
#####################################

# depends -----------------------------------------------------------------

options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
suppressMessages(suppressWarnings(if (!require('getopt')) install.packages('getopt')));
suppressMessages(suppressWarnings(if (!require('logr')) install.packages('logr')));
# parameters --------------------------------------------------------------
tmp <- file.path("Step1.log")

lf <- log_open(tmp)

command=matrix(c(
  'help', 'h', 0, 'logic', 'help information',
  'datExpr', 'd', 1, 'character', 'inputfile: expression matrix, geneID in row, sample names in column',
  'datatype', 't', 1, 'character', 'Data type: count or FPKM. count means readcount, FPKM represent normalized count,eg: FPKM TPM CPM RPKM, Default: count',
  'method', 'm', 1, 'character', 'Data transformat method: For count: varianceStabilizingTransformation or cpm; For FPKM, rawFPKM or logFPKM, Default: varianceStabilizingTransformation',
  'RcCutoff', 'r', 1, 'double', 'Noise cutoff: count/normalized count lower than ... was thought to be noise. Default: count=>10, Normalized count=>1',
  'samplePerc', 's', 1, 'double', 'Sample percentage: parameter for noise remove. xx percent of all samples have readcount/FPKM > cutoff. Default: 0.3',
  'GeneNum', 'g', 1, 'double', 'Gene number for WGCNA: how many gene you want retained for WGCNA after noise remove.',
  'cutmethod', 'a', 1, 'character', 'MAD or SVR, Default: MAD',
  'rscut', 'b', 1, 'double', 'sft Power cutoff: Power cutoff',
  'threads', 'c', 1, 'double', 'Threads, how many threads you want to use,Default = 8'
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

if (is.null(args$threads)){
  args$threads = 8
}
threads = args$threads


## default value
if (is.null(args$datExpr)){
  q(status = 1)
  log_print(
    "!<== ERROR: no input file. ==>"
  )
}

datExpr <- read.table(
  args$datExpr,header = T,sep = "\t"
)

ngenes = nrow(datExpr)

if (is.null(args$datatype)){
  args$datatype = 'count'
}

datatype = args$datatype

if (is.null(args$method)){
  args$method = "varianceStabilizingTransformation"
}

method = args$method

if (is.null(args$RcCutoff)){
  args$RcCutoff = 10
}

RcCutoff = args$RcCutoff

if (is.null(args$samplePerc)){
  args$samplePerc = 0.3
}

samplePerc = args$samplePerc

if (is.null(args$GeneNum)){
  args$GeneNum = 8000
}
GeneNum = args$GeneNum
GeneNumCut = 1 - (GeneNum/ngenes)

if (is.null(args$cutmethod)){
  args$cutmethod = "MAD"
}
cutmethod = args$cutmethod
if (is.null(args$rscut)){
  args$rscut = 0.8
}

rscut = args$rscut
# run ---------------------------------------------------------------------
suppressMessages(suppressWarnings(if (!requireNamespace("BiocManager", quietly = TRUE))
				    install.packages("BiocManager")))
suppressMessages(suppressWarnings(if (!require('devtools')) install.packages('devtools')));
suppressMessages(suppressWarnings(if (!require('WGCNA')) BiocManager::install('WGCNA',update = FALSE)));
suppressMessages(suppressWarnings(if (!require('tidyverse')) install.packages('tidyverse')));
suppressMessages(suppressWarnings(if (!require('ShinyWGCNA')) devtools::install_github("ShawnWx2019/WGCNAShinyFun",ref = "master")));
suppressMessages(suppressWarnings(if (!require('DESeq2')) BiocManager::install('DESeq2',update = FALSE)));
suppressMessages(suppressWarnings(if (!require('ggprism')) BiocManager::install('ggprism',update = FALSE)));
suppressMessages(suppressWarnings(if (!require('patchwork')) BiocManager::install('patchwork',update = FALSE)));
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
  "<== Noise remove ==>"
)

tryCatch({
  datExpr <-
    getdatExpr(
      rawdata = datExpr,
      RcCutoff = RcCutoff,
      samplePerc = samplePerc,
      datatype = datatype,
      method = method
    )
  log_print(
    "Step1. Noise remove is Success!"
  )
},error = function(e) {
  log_print(
    "Step1. Noise remove is Failed!"
  )
  q(status = 1)
})

log_print(
  "<== Step2 filter ==>"
)


tryCatch(
  {
    datExpr <-
      getdatExpr2(
        datExpr = datExpr,
        GeneNumCut = GeneNumCut,
        cutmethod = cutmethod
      )
    log_print(
      "Step2. Gene filtering is Success!"
    )
  },error = function(e) {
    log_print(
      "Step2. Gene filtering is Failed!"
    )
    q(status = 1)
  }
)


tryCatch(
  {
    step1_sft = getpower(
      datExpr = datExpr,
      rscut = rscut
    )
    log_print(
     paste0( "Step3. SFT detective is Success, The recommanded power is: ", step1_sft$power)
    ) 
    log_print(
	      step1_sft$sft
    )
  },error = function(e) {
    log_print(
      "Step3. SFT detective is Failed!"
    )

  }
)
params = data.frame(
	param = c("datatype","method","RcCutoff","samplePerc","GeneNum","cutmethod","rscut","threads"),
	value = c(datatype,method,RcCutoff,samplePerc,GeneNum,cutmethod,rscut,threads)
)

log_print("Paramters used in this step: ")
log_print(params)

save.image("step1.RData")

log_print(
  "< == All work finish ==>"
)

log_close()

writeLines(readLines(lf))

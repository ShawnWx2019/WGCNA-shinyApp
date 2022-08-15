#####################################
# 	Desc: WGCNA for server pipeline
# 	Date: Aug 11, 2022
#   File: Step2. module-trait
#	  Author: Shawn Wang
#####################################

# depends -----------------------------------------------------------------
suppressMessages(suppressWarnings(library(getopt)));
suppressMessages(suppressWarnings(library(logr)));
load("step2.RData")
options(scipen = 6)
# parameters --------------------------------------------------------------
tmp <- file.path("Step3.log")

lf <- log_open(tmp)

command=matrix(c(
  'help', 'h', 0, 'logic', 'help information',
  'traitData', 't' , 1, 'character', 'Trait data, Sample ids in row and traits in column.'
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
if (is.null(args$traitData)){
  q(status = T)
}


# run ---------------------------------------------------------------------
suppressMessages(suppressWarnings(library(WGCNA)));
suppressMessages(suppressWarnings(library(ShinyWGCNA)));
suppressMessages(suppressWarnings(library(tidyverse)));
suppressMessages(suppressWarnings(library(ggprism)));
suppressMessages(suppressWarnings(library(patchwork)));
suppressMessages(suppressWarnings(if (!require('ComplexHeatmap')) BiocManager::install('ComplexHeatmap',update = FALSE)));
suppressMessages(suppressWarnings(if (!require('circlize')) BiocManager::install('circlize',update = FALSE)));
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

fpath <- args$traitData

if(str_detect(fpath,".csv")) {
  phen <- read.csv(fpath)
} else if(str_detect(fpath,".xlsx")) {
  phen <- readxl::read_xlsx(fpath)
} else {
  phen <- read.table(fpath,header = T,sep = "\t")
}
colnames(phen)[1] = "sample_id"
phen <- left_join(
  data.frame(
    sample_id = rownames(datExpr)
  ),phen,"sample_id"
) %>% 
  column_to_rownames("sample_id")


tryCatch({
  
  step3_mdule_trait <-
    getMt(
      phenotype = phen,
      nSamples = nrow(datExpr),
      moduleColors = step2_network$moduleColors,
      datExpr = datExpr
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

mod_color = gsub(pattern = "^..",replacement = "",rownames(step3_mdule_trait$modTraitCor))
mod_color_anno = setNames(mod_color,rownames(step3_mdule_trait$modTraitCor))


left_anno = rowAnnotation(
  Module = rownames(step3_mdule_trait$modTraitCor),
  col = list(
    Module = mod_color_anno
  ),
  show_legend = F,
  show_annotation_name = F
)

m_t_ht <- 
  Heatmap(
    matrix = step3_mdule_trait$modTraitCor,
    cluster_rows = F, cluster_columns = F,
    left_annotation = ,
    cell_fun = function(j,i,x,y,width,height,fill) {
      grid.text(sprintf(step3_mdule_trait$textMatrix[i,j]),x,y,gp = gpar(fontsize = 12))
    },
    row_names_side = "left",
    column_names_rot = 45,
    heatmap_legend_param = list(
      at = c(-1,-0.5,0,0.5, 1),
      labels = c("-1","-0.5", "0","0.5", "1"),
      title = "",
      legend_height = unit(9, "cm"),
      title_position = "lefttop-rot"
    ),
    rect_gp = gpar(col = "black", lwd = 1.2),
    column_title = "Module-trait relationships",
    column_title_gp = gpar(fontsize = 15, fontface = "bold"),
    col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
  )

pdf("module-trait.pdf",width = 18,height = 12)
draw(m_t_ht)
dev.off
params = data.frame(
  param = c("datatype","method","RcCutoff","samplePerc","GeneNum","cutmethod","rscut","threads","power","minModuleSize","mergeCutHeight","maxBlocksize"),
  value = c(datatype,method,RcCutoff,samplePerc,GeneNum,cutmethod,rscut,threads,power,minModuleSize,mergeCutHeight,maxBlocksize)
)

log_print("Paramters used in this step: \n")
log_print(params)

save.image("step3.RData")

log_print(
  "< == All work finish ==>"
)

log_close()

writeLines(readLines(lf))

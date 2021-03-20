# fuction ===============================
getdatExpr = function(rawdata,RcCutoff,samplePerc,datatype,method){
  ## Avoid numeric id 
  rawdata <- data.frame(row.names = as.character(rawdata[,1]),
                        rawdata[,-1])
  countvarTran = function(rawcount,RcCutoff,samplePerc) {
    samnum <- ncol(rawcount)
    casenum = ceiling(samnum/2)
    controlnum = samnum - casenum
    condition <- factor(c(rep("case",casenum),rep("control",controlnum)),levels = c("case","control"))
    colData <- data.frame(row.names = colnames(rawcount), condition)
    ## remove background noise
    x <- rawcount[apply(rawcount,1,function(x) sum(x > RcCutoff) > (samplePerc*ncol(rawcount))),]
    
    ## readcount standardization by DESeq2
    dds <- DESeqDataSetFromMatrix(x, colData, design = ~ condition)
    dds <- DESeq(dds)
    vsd <- assay(varianceStabilizingTransformation(dds))
    dx = data.frame(vsd)
    return(dx)
  }
  ## func2 count 2 cpm
  countCPM = function(rawcount, RcCutoff,samplePerc) {
    x <- rawcount[apply(rawcount,1,function(x) sum(x > RcCutoff) > (samplePerc*ncol(rawcount))),]
    dx  <- log10(edgeR::cpm(x)+1)
    return(dx)
  }
  ## raw fpkm filter
  fpkmfilter = function(rawcount, RcCutoff,samplePerc) {
    x <- rawcount[apply(rawcount,1,function(x) sum(x > RcCutoff) > (samplePerc*ncol(rawcount))),]
    dx  <- x
    return(dx)
  }
  
  ## log fpkm filter
  lgfpkmfilter = function(rawcount, RcCutoff,samplePerc) {
    x <- rawcount[apply(rawcount,1,function(x) sum(x > RcCutoff) > (samplePerc*ncol(rawcount))),]
    dx  <- log10(x+1)
    return(dx)
  }
  ## dx final
  if (
    datatype == "count" & method == "varianceStabilizingTransformation"
  ) {
    dx = countvarTran(rawcount = rawdata,RcCutoff = RcCutoff, samplePerc = samplePerc)
  } else if (
    datatype == "count" & method == "lgcpm"
  ) {
    dx = countCPM(rawcount = rawdata,RcCutoff = RcCutoff, samplePerc = samplePerc)
  } else if (
    datatype == "FPKM" & method == "rawFPKM"
  ) {
    dx = fpkmfilter(rawcount = rawdata,RcCutoff = RcCutoff, samplePerc = samplePerc)
  } else if (
    datatype == "FPKM" & method == "lgFPKM"
  ) {
    dx = lgfpkmfilter(rawcount = rawdata,RcCutoff = RcCutoff, samplePerc = samplePerc)
  }
  return(dx)
}
getdatExpr2 = function(datExpr,GeneNumCut,cutmethod){
  datExpr <- datExpr
  type = "unsigned"
  corType = "pearson"
  corFnc = cor
  maxPOutliers = ifelse(corType=="pearson",1,0.05)
  robustY = ifelse(corType=="pearson",T,F)
  m.mad <- apply(datExpr,1,mad)
  m.var <- apply(datExpr,1,var)
  if(GeneNumCut == 0){
    datExprVar = datExpr
  } else if(cutmethod == "MAD"){
    datExprVar <- datExpr[which(m.mad > 
                                  max(quantile(m.mad, probs=seq(0, 1, GeneNumCut))[2],0.01)),]
  } else if(cutmethod == "Var"){
    datExprVar <- datExpr[which(m.var > 
                                  max(quantile(m.var, probs=seq(0, 1, GeneNumCut))[2],0.01)),]
  } else {
    print("Error: wrong filter method")
  }
  dim(datExprVar)
  datExpr <- as.data.frame(t(datExprVar))
  gsg = goodSamplesGenes(datExpr, verbose = 3)
  if (!gsg$allOK){
    # Optionally, print the gene and sample names that were removed:
    if (sum(!gsg$goodGenes)>0) 
      printFlush(paste("Removing genes:", 
                       paste(names(datExpr)[!gsg$goodGenes], collapse = ",")));
    if (sum(!gsg$goodSamples)>0) 
      printFlush(paste("Removing samples:", 
                       paste(rownames(datExpr)[!gsg$goodSamples], collapse = ",")));
    # Remove the offending genes and samples from the data:
    datExpr = datExpr[gsg$goodSamples, gsg$goodGenes]
  }
  return(datExpr)
}
getsampleTree = function(datExpr,layout = "circular"){
  ## sample cluster based on expression values
  nGenes = ncol(datExpr)
  nSamples = nrow(datExpr)
  assign("nGenes",value = nGenes, envir = globalenv())
  assign("nSamples",value = nSamples, envir = globalenv())
  sampleTree = hclust(dist(datExpr), method = "average")
  treenew= as.phylo(sampleTree)
  anno = data.frame(lab = treenew$tip.label,
                    group = factor(gsub(".$","",treenew$tip.label),levels = unique(gsub(".$","",treenew$tip.label))))
  p =ggtree(tr = treenew,layout = layout) 
  p2 = p %<+% anno +geom_tiplab(aes(color=group))+theme(legend.position = "")
  out = list(tree = treenew,
             plot = p2)
  return(out)
}
getpower = function(datExpr,rscut){
  type = 'unsigned'
  powers = c(c(1:10), seq(from = 12, to=30, by=2))
  sft = pickSoftThreshold(datExpr, powerVector=powers, RsquaredCut = rscut,
                          networkType=type, verbose=5)
  fitIn = sft$fitIndices
  min.sft = min(-sign(sft$fitIndices[,3])*sft$fitIndices[,2])
  fitIn %>% 
    mutate(r2 = -sign(slope)*SFT.R.sq) %>% 
    ggplot(.,aes(x = Power,y = r2,color = "red"))+
    geom_text(aes(label = Power))+
    theme_bw()+
    ylim(min.sft,1)+
    geom_hline(aes(yintercept = rscut),color = "green",linetype = "dashed")+
    geom_hline(aes(yintercept = 0.8),color = "blue",linetype = "dashed")+
    geom_hline(aes(yintercept = 0.9),color = "red",linetype = "solid")+
    theme_prism(border = T)+
    xlab("Soft Threshold (Power)")+
    ylab(expression("Scale Free Topology Model Fit,signed R"^2))+
    ggtitle("Scale independence")+
    theme(legend.position = "")->p1
  fitIn %>% 
    ggplot(.,aes(x = Power,y = mean.k.,color = "red"))+
    geom_text(aes(label = Power))+
    theme_bw()+
    theme_prism(border = T)+
    xlab("Soft Threshold (Power)")+
    ylab("Mean Connectivity")+
    ggtitle("Mean Connectivity")+
    theme(legend.position = "")->p2
  p.all = p1+p2
  power = sft$powerEstimate
  ## exper power
  if (is.na(power)){
    power = ifelse(nSamples<20, ifelse(type == "unsigned", 9, 18),
                   ifelse(nSamples<30, ifelse(type == "unsigned", 8, 16),
                          ifelse(nSamples<40, ifelse(type == "unsigned", 7, 14),
                                 ifelse(type == "unsigned", 6, 12))
                   )
    )
  }
  out = list(sft = fitIn,
             plot = p.all,
             power = power)
  
  return(out)
}
ggscaleFreePlot = function (connectivity, nBreaks = 10, truncated = FALSE, removeFirst = FALSE, main = "", ...) {
  k = connectivity
  discretized.k = cut(k, nBreaks)
  dk = tapply(k, discretized.k, mean)
  p.dk = as.vector(tapply(k, discretized.k, length)/length(k))
  breaks1 = seq(from = min(k), to = max(k), length = nBreaks + 
                  1)
  hist1 = suppressWarnings(hist(k, breaks = breaks1, equidist = FALSE, 
                                plot = FALSE, right = TRUE, ...))
  dk2 = hist1$mids
  dk = ifelse(is.na(dk), dk2, dk)
  dk = ifelse(dk == 0, dk2, dk)
  p.dk = ifelse(is.na(p.dk), 0, p.dk)
  log.dk = as.vector(log10(dk))
  if (removeFirst) {
    p.dk = p.dk[-1]
    log.dk = log.dk[-1]
  }
  log.p.dk = as.numeric(log10(p.dk + 1e-09))
  lm1 = lm(log.p.dk ~ log.dk)
  if (truncated == TRUE) {
    lm2 = lm(log.p.dk ~ log.dk + I(10^log.dk))
    OUTPUT = data.frame(scaleFreeRsquared = round(summary(lm1)$adj.r.squared, 
                                                  2), slope = round(lm1$coefficients[[2]], 2), TruncatedRsquared = round(summary(lm2)$adj.r.squared, 
                                                                                                                         2))
    printFlush("the red line corresponds to the truncated exponential fit")
    title = paste(main, " scale free R^2=", as.character(round(summary(lm1)$adj.r.squared, 
                                                               2)), ", slope=", round(lm1$coefficients[[2]], 2), 
                  ", trunc.R^2=", as.character(round(summary(lm2)$adj.r.squared, 
                                                     2)))
  }
  else {
    title = paste(main, " scale R^2=", as.character(round(summary(lm1)$adj.r.squared, 
                                                          2)), ", slope=", round(lm1$coefficients[[2]], 2))
    OUTPUT = data.frame(scaleFreeRsquared = round(summary(lm1)$adj.r.squared, 
                                                  2), slope = round(lm1$coefficients[[2]], 2))
  }
  pdata = data.frame(x = log.dk,
                     y = log.p.dk)
  lm = if(truncated == F){
    lm1
  } else {
    lm2
  }
  
  p = ggplot(pdata ,aes(x = x,y = y))+
    geom_point(shape  = 1,size = 5)+
    geom_abline(intercept = as.numeric(lm$coefficients[1]),slope =as.numeric(lm$coefficients[2]) )+
    xlab("log10(k)")+
    ylab("log10(p(k))")+
    ggtitle(title)+theme_prism(border = T)
  
}
powertest = function(power.test,datExpr,nGenes){
  ## reference: https://www.jianshu.com/p/25905a905086 jiashu Author: sxzt
  softPower <- power.test
  adjacency = adjacency(datExpr = datExpr, power = softPower)
  # convert adj 2 TOM
  TOM = TOMsimilarity(adjacency)
  # disTOM
  dissTOM = 1-TOM
  hierTOM = hclust(as.dist(dissTOM),method="average")
  ADJ1_cor <- abs(WGCNA::cor( datExpr,use = "p" ))^softPower
  if(nGenes< 5000){
    k <- as.vector(apply(ADJ1_cor,2,sum,na.rm=T))
  } else {
    k <- softConnectivity(datE=datExpr,power=softPower) 
  }
  k.df = data.frame(value = k)
  p1 = ggplot(k.df,mapping = aes(x = value)) + 
    geom_histogram(color = "black",fill = "white")+
    theme_prism(border = T)+
    xlab("k")+
    ylab("Frequence")+
    ggtitle(paste("soft connectivity (power = ",power.test,")"))+
    theme(legend.position = "")
  p2 = ggscaleFreePlot(k,main = "Check Scale free topology\n")
  p = p1+p2
  return(p)
}
getnetwork = function(datExpr,power,minModuleSize,mergeCutHeight){
  cor <- WGCNA::cor
  net = blockwiseModules(datExpr, power = power, maxBlockSize = nGenes,
                         TOMType = "unsigned", minModuleSize = minModuleSize,
                         reassignThreshold = 0, mergeCutHeight = mergeCutHeight,
                         numericLabels = TRUE, pamRespectsDendro = FALSE,
                         saveTOMs=F, corType = 'pearson', 
                         maxPOutliers=1,
                         verbose = 3)
  moduleLabels = net$colors
  moduleColors = labels2colors(moduleLabels)
  MEs = net$MEs
  MEs_col = MEs
  colnames(MEs_col) = paste0("ME", labels2colors(
    as.numeric(str_replace_all(colnames(MEs),"ME",""))))
  MEs_col = orderMEs(MEs_col)
  x = table(moduleColors)
  Gene2module <- data.frame(GID = colnames(datExpr),
                            Module = moduleColors)
  out = list(net = net,
             moduleLabels =moduleLabels,
             moduleColors = moduleColors,
             MEs_col = MEs_col,
             MEs = MEs,
             Gene2module = Gene2module
  )
  return(out)
}
getMt = function(phenotype,MEs_col,nSamples,moduleColors,datExpr){
  traitData <- phenotype
  # if (corType=="pearson") {
  #   modTraitCor = cor(MEs_col, traitData, use = "p")
  #   modTraitP = corPvalueStudent(modTraitCor, nSamples)
  # } else {
  #   modTraitCorP = bicorAndPvalue(MEs_col, traitData, robustY=robustY)
  #   modTraitCor = modTraitCorP$bicor
  #   modTraitP   = modTraitCorP$p
  # }
  MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
  MEs = orderMEs(MEs0)
  modTraitCor = cor(MEs, traitData , use = "p")
  modTraitP = corPvalueStudent(modTraitCor, nSamples)
  
  textMatrix = paste(signif(modTraitCor, 2), 
                     "\n(", signif(modTraitP, 1), ")", sep = "")
  
  out = list(modTraitCor = modTraitCor,
             modTraitP = modTraitP,
             textMatrix = textMatrix)
  return(out)
  
}
getKME = function(datExpr,moduleColors,MEs_col){
  connet=abs(cor(datExpr,use="p"))^6
  Alldegrees1=intramodularConnectivity(connet, moduleColors)
  datKME=signedKME(datExpr, MEs_col, outputColumnName="MM.")
  return(datKME)
}
getMM = function(datExpr,MEs_col,nSamples,corType){
  if (corType=="pearson") {
    geneModuleMembership = as.data.frame(cor(datExpr, MEs_col, use = "p"))
    MMPvalue = as.data.frame(corPvalueStudent(
      as.matrix(geneModuleMembership), nSamples))
  } else {
    geneModuleMembershipA = bicorAndPvalue(datExpr, MEs_col, robustY=robustY)
    geneModuleMembership = geneModuleMembershipA$bicor
    MMPvalue   = geneModuleMembershipA$p
  }
  out = list(MM= geneModuleMembership,
             MMP = MMPvalue)
}
getverboseplot = function(datExpr,module,pheno,traitData,moduleColors,geneModuleMembership,MEs,nSamples){
  modNames = substring(names(MEs), 3)
  tmp1 = traitData %>% 
    select(pheno);
#  print(dim(tmp1))
#  print(dim(datExpr))
  geneTraitSignificance = as.data.frame(cor(datExpr, tmp1, use = "p"));
#  print(dim(geneTraitSignificance))
  GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
  names(geneTraitSignificance) = paste("GS.", names(tmp1), sep="");
#  print(dim(geneTraitSignificance))
  names(GSPvalue) = paste("p.GS.", names(tmp1), sep="")
  module = module
  column = match(module, modNames);
  moduleGenes = moduleColors==module;
  plt.df = data.frame(x = abs(geneModuleMembership[moduleGenes, column]),
                      y = abs(geneTraitSignificance[moduleGenes, 1]))
  #formula <- y ~ x
  corExpr = parse(text = paste("cor", "(plt.df$x, plt.df$y ", prepComma("use = 'p'"), 
                               ")"))
  cor = signif(eval(corExpr), 2)
  if (is.finite(cor)) {
    if (abs(cor) < 0.00001){ cor = 0 } 
  }
  corp = signif(corPvalueStudent(cor, sum(is.finite(plt.df$x) & is.finite(plt.df$y))), 2)
  if (is.finite(corp) && corp < 10^(-200)) {
    corp = "<1e-200"
  } else {corp = paste("=", corp, sep = "")}
  
  ggplot(plt.df,aes(x = x, y = y))+
    geom_point(shape = 21,color = "black",fill = module,size = 3)+
    geom_smooth(method = 'lm', formula = y ~ x)+
    xlab(paste("Module Membership in", module, "module"))+
    ylab(paste("Gene significance for ",pheno))+
    ggtitle(paste("Module membership vs. gene significance\n","cor=",cor," ,p",corp,sep = ""))+
    theme_bw()+
    theme_prism(border = T)

  
  # verboseScatterplot(x = abs(geneModuleMembership[moduleGenes, column]),abline = T,
  #                    y = abs(geneTraitSignificance[moduleGenes, 1]),
  #                    xlab = paste("Module Membership in", module, "module"),
  #                    ylab = paste("Gene significance for ",pheno),
  #                    main = paste("Module membership vs. gene significance\n"),
  #                    cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
}
moduleheatmap = function(datExpr,MEs,which.module,moduleColors){
  ME=MEs[,paste("ME",which.module,sep = "")]
  ## heatmap data
  h <- t(scale(datExpr[,moduleColors == which.module]))
  x = reshape2::melt(data = h,id.var = rownames(h))
  ## barplot data
  b <- data.frame(Sample = factor(rownames(datExpr),levels = rownames(datExpr)) ,
                  Eigengene = ME)
  p1 = ggplot(x, aes(Var2, Var1)) + 
    geom_tile(aes(fill = value),colour = alpha("white",0)) + 
    scale_fill_gradient2(low = "green",mid = "black",high = "red")+
    theme_classic()+
    theme(legend.position = "none",
          axis.title = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_text(angle = 90,face = "bold",vjust = 0.5,size = 15),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.x = element_blank())
  ## barplot
  p2 = ggplot(b,mapping  = aes(x = Sample,y = Eigengene)) + geom_bar(stat='identity',fill = which.module)+
    theme_classic()+
    theme(axis.title = element_blank(),
          axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())
  p3 = p1/p2
  return(p3)
}
## hubgenen
hubgenes = function(datExpr,mdl,power,trt,KME,GS.cut,kME.cut,datTrait){
  KME.clean = KME
  colnames(KME.clean) = gsub("MM.","",colnames(KME))
  hub1 = chooseTopHubInEachModule(datExpr = datExpr,colorh = mdl,power = power,type = "unsigned")
  ## GS calculate
  which.module = mdl
  which.trait = datTrait %>% 
    select(trt); 
  GS1 = as.numeric(cor(which.trait,datExpr, use="p"))
  GS.df = data.frame(GeneID = colnames(datExpr),
                     GS = GS1)
  GS2=abs(GS1)
  ## KME
  kME.mod = KME.clean %>% 
    select(mdl)
  kME.mod.df = data.frame(GeneID = rownames(kME.mod),
                          kME = kME.mod[,1])
  ## judge
  hub2 = abs(kME.mod) > kME.cut & GS2 > kME.cut
  hub2.num = table(hub2)
  print(paste0("kME and GS method selected out " ,as.numeric(hub2.num[2])," hub genes."))
  colnames(hub2) = "judge"
  hub3 <- hub2 %>%as.data.frame() %>%  filter(judge == 'TRUE') %>% 
    rownames_to_column(var = "GeneID") %>% 
    left_join(.,GS.df,by = "GeneID") %>% 
    left_join(.,kME.mod.df,by = "GeneID")
  hub1 = data.frame(GeneID = as.character(hub1)) %>% 
    left_join(.,GS.df,by = "GeneID") %>% 
    left_join(.,kME.mod.df,by = "GeneID")


  out = list(hub1 = hub1,
             hub3 = hub3)
}

## 
TOM = TOMsimilarityFromExpr(datExpr, power = sft$powerEstimate); 
# Select module
module = "turquoise";
# Select module probes
probes = colnames(datExpr) ## 我们例子里面的probe就是基因名
inModule = (moduleColors==module);
modProbes = probes[inModule]; 
## 也是提取指定模块的基因名
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)

cyt = exportNetworkToCytoscape(
  modTOM,
  edgeFile = paste("CytoscapeInput-edges-", paste(module, collapse="-"), ".txt", sep=""),
  nodeFile = paste("CytoscapeInput-nodes-", paste(module, collapse="-"), ".txt", sep=""),
  weighted = TRUE,
  threshold = 0.02,
  nodeNames = modProbes, 
  nodeAttr = moduleColors[inModule]
);

# run ===========================
x = read.delim("~/03.project/01.Heterosis/03.process/01.S9108/04.WGCNA/fiber.fpkm.txt",header = T,sep = "\t",stringsAsFactors = F)
y = getdatExpr(rawdata  = x,RcCutoff = 1,samplePerc = 0.3,datatype = "FPKM",method = "rawFPKM" )
### filter step2 ===================================================

RemainGeneNum = 1000
GNC = 1-RemainGeneNum/nrow(y)
z = getdatExpr2(datExpr = y,GeneNumCut = GNC,cutmethod = "MAD")
### Sample tree===============================

getsampleTree(datExpr = z)
##Part2 sft=========================

### getpower==============


tmp1 = getpower(datExpr = z,rscut = 0.9)
## testpower============================
### change scaleFreePlot into gg format
tmp1$plot
tmp1$power
datExpr = z
powertest(power.test = 14,datExpr = z,nGenes = ncol(z))

dd= getnetwork(datExpr = z,power = tmp1$power,minModuleSize = 30,mergeCutHeight = 0.25)
#
plotDendroAndColors(dd$net$dendrograms[[1]], dd$moduleColors[dd$net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

trait = read.delim("~/03.project/01.Heterosis/03.process/01.S9108/04.WGCNA/fiber.trait.txt",
                       header = T,sep = "\t")

if (ncol(trait) == 2) {
  x <- trait
  Tcol = as.character(unique(x[,2]))
  b <- list()
  for (i in 1:length(Tcol)) {
    b[[i]] = data.frame(row.names = x[,1],
                        levels = ifelse(x[,2] == Tcol[i],1,0))
  }
  c <- bind_cols(b)
  c <- data.frame(row.names = x$name,
                  c)
  colnames(c) = Tcol
  rownames(c) = trait[,1]
  pheTmp <- c
} else {
  pheTmp = data.frame(row.names = trait[,1],
                      trait[,-1])
}
exp.ds = getMt(phenotype = pheTmp,MEs_col = dd$MEs_col,nSamples = nrow(z),
               moduleColors = dd$moduleColors,datExpr = z)
library(tidyverse)
## getKME


labeledHeatmap(Matrix = exp.ds$modTraitCor, xLabels = colnames(pheTmp), 
               yLabels = colnames(dd$MEs_col), 
               cex.lab = 0.7, xLabelsAngle = 45, xLabelsAdj = 1,
               ySymbols = substr(colnames(dd$MEs_col),3,1000), colorLabels = FALSE, 
               colors = colorRampPalette(c("orange","white","purple"))(100), 
               textMatrix = exp.ds$textMatrix, setStdMargins = FALSE, 
               cex.text = 0.6, zlim = c(-1,1),
               main = paste("Module-trait relationships"))
KME = getKME(datExpr = z,moduleColors = dd$moduleColors,MEs_col = dd$MEs_col)

mm.out = getMM(datExpr = z,MEs_col = dd$MEs_col,nSamples = ncol(z),corType = "pearson")


head(KME)
## 导出hubgenen
KME.clean = KME
colnames(KME.clean) = gsub("MM.","",colnames(KME))
hub.gene = hubgenes(datExpr = z,mdl = "green",power = 22,trt = "Z41S",KME = KME,GS.cut = 0.5,kME.cut = 0.5,datTrait = pheTmp)
hub.gene$hub3
save(getdatExpr,getdatExpr2,
     getsampleTree,getpower,
     ggscaleFreePlot,powertest,
     getnetwork,getMt,getKME,getMM,
     getverboseplot,moduleheatmap,hubgenes,
     file = "~/02.MyScript/OneStepWGCNA/03.shinyAPP/functions.Rdata")



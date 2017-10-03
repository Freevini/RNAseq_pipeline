args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

## This script performs differential expression analysis with edgeR, based on 
## abundance estimates from Salmon. It supports testing one or more contrasts. 
## To be compatible with downstream analysis scripts, the output has to be
## either:
## - a data frame with results (e.g., returned by edgeR's topTags(...)$table).
## There must be one column named "gene" that gives the gene ID.
## - a named list of such data frames, one for each contrast (recommended).
## To run the script, modify at least the definition of the design matrix and
## the contrasts of interest.

##variables to do it local
#tx2gene <- "/home/bioinf/bioinf_data/43_sovi/Projects/Metagenomics_pipeline/Model_data/e_coli/rnaseqworkflow-master/reference/Escherichia_coli_o157_h7_str_sakai.ASM886v1.cdna.ncrna_tx2gene.rds"
#salmondir <- "/home/bioinf/bioinf_data/43_sovi/Projects/Metagenomics_pipeline/Model_data/e_coli/rnaseqworkflow-master/salmon"
#metafile <- "/home/bioinf/bioinf_data/43_sovi/Projects/Metagenomics_pipeline/Model_data/e_coli/rnaseqworkflow-master/metadata.txt"
# outrds <- "/home/bioinf/bioinf_data/43_sovi/Projects/Metagenomics_pipeline/Model_data/e_coli/rnaseqworkflow-master/output/edgeR_dge.rds"
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tximport))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library("vegan"))

# source("/home/bioinf/bioinf_data/43_sovi/Projects/Metagenomics_pipeline/Model_data/e_coli/05_trial_05/scripts/mean_variance_script.R")
source("scripts/mean_variance_script.R")


print(tx2gene)
print(salmondir)
print(metafile)
print(outrds)

## List Salmon directories
salmondirs <- list.files(salmondir, full.names = TRUE)
salmonfiles <- paste0(salmondirs, "/quant.sf")
names(salmonfiles) <- basename(salmondirs)
(salmonfiles <- salmonfiles[file.exists(salmonfiles)])

## Read transcript-to-gene mapping
tx2gene <- readRDS(tx2gene)

## Read Salmon abundances
txi <- tximport(files = salmonfiles, type = "salmon", txOut = FALSE, 
                tx2gene = tx2gene[, c("tx", "gene")],dropInfReps=TRUE)

## Read metadata and reorder in the same order as the count matrix
metadata <- read.delim(metafile, header = TRUE, as.is = TRUE, sep = "\t")
stopifnot(all(metadata$ID %in% colnames(txi$counts)))
stopifnot(all(colnames(txi$counts) %in% metadata$ID))
metadata <- metadata[match(colnames(txi$counts), metadata$ID), ]

## Create DGEList and include average transcript length offsets
cts <- txi$counts
normMat <- txi$length
normMat <- normMat/exp(rowMeans(log(normMat)))
o <- log(calcNormFactors(cts/normMat)) + log(colSums(cts/normMat))
dge0 <- DGEList(cts)
dge0$offset <- t(t(log(normMat)) + o)
dge0 <- calcNormFactors(dge0)

## Define design. ************** MODIFY ************** 
stopifnot(all(colnames(dge0) == metadata$ID))
(des <- model.matrix(~ time, data = metadata))

## Filter out genes with average CPM below 1
print(dim(dge0))
cpms <- cpm(dge0)
dge <- dge0[apply(cpms, 1, mean) > 1, ]
dge <- calcNormFactors(dge)
print(dim(dge))

## Add gene annotation
annot <- tx2gene %>% dplyr::select(-tx, -tx_biotype, -start, -end) %>% distinct()
if (any(duplicated(annot$gene))) {
  stop(paste0("The following genes are represented by multiple rows in the ", 
              "gene annotation: ", 
              paste(annot$gene[duplicated(annot$gene)], collapse = ",")))
}
annot <- annot[match(rownames(dge), annot$gene), ]
rownames(annot) <- annot$gene
dge$genes <- annot

dge_cp <- dge
## Estimate dispersion and fit model
dge <- estimateDisp(dge, design = des) ## calculate all common,tagwise and trended dispersion
theta <- 1/dge$common.dispersion
qlfit <- glmQLFit(dge, design = des)

##----
###plot mean-variance (Oberg et al. 2012)
##----
jpeg(filename = "diagnosticPlots/Mean-variance-plot.jpeg")
plotMeanVariance(dge0$counts,metadata$time)
dev.off()

##----
##plot goodness of fit (GOF)
##----

##global 
gcom_true <- gof(qlfit)
zcom_true <- zscoreGamma(gcom_true$gof.statistics,shape=gcom_true$df/2,scale=2)
pcom_ture=zcom_true[is.finite(zcom_true) ]


y_1<-estimateGLMCommonDisp(dge_cp, design = des) 
fitcommon <- glmQLFit(y_1,design = des) 
y_2 <- estimateGLMTrendedDisp(dge_cp, design = des) 
fittrended <- glmQLFit(y_2,design = des) 
y_3 <- estimateGLMTagwiseDisp(y_1, design = des) 
fittagwise <- glmQLFit(y_3,design = des)

gcom <- gof(fitcommon) 
gtren <- gof(fittrended) 
gtag <- gof(fittagwise)


zcom <- zscoreGamma(gcom$gof.statistics,shape=gcom$df/2,scale=2) 
ztren <- zscoreGamma(gtren$gof.statistics,shape=gtren$df/2,scale=2) 
ztag <- zscoreGamma(gtag$gof.statistics,shape=gtag$df/2,scale=2)


pcom=zcom[is.finite(zcom) ] 
ptren=ztren[is.finite(ztren) ] 
ptag=ztag[is.finite(ztag) ]

##jpeg(filename = "/home/bioinf/bioinf_data/43_sovi/Projects/Metagenomics_pipeline/Model_data/e_coli/06_trial_06/diagnosticPlots/Q-Q-plots.jpeg")
jpeg(filename = "diagnosticPlots/Q-Q-plots.jpeg")

par(mfrow=c(2,2))

qqnorm(pcom_ture,main = "weighted_likelihood-empirical_Bayes-Q-Q-plot")
qqline(pcom_ture)

qqnorm(pcom,main = "common-Q-Q-plot")
qqline(pcom) 

qqnorm(ptren,main = "trended-Q-Q-plot") 
qqline(ptren) 

qqnorm(ptag,main = "tagwise-Q-Q-plot") 
qqline(ptag)
dev.off()

##----
##Rpkm
##----
genLength <- (as.integer(tx2gene$end)-as.integer(tx2gene$start))
Rpkm <- rpkm(dge,gene.length = genLength)

##----
##average count per gene
##----
jpeg(filename = "diagnosticPlots/count_histogram.jpeg")
# jpeg(filename = "diagnosticPlots/count_histogram.jpeg")
par(mfrow=c(1,1))
a <- 100 ##number of breaks
hist(log10(rowMeans(dge$counts)),main="Frequency histogram of average counts per gene",
     xlab="log10(Average counts per gene)",
     ylab="Frequency",
     breaks = seq(from=min(log10(rowMeans(dge$counts))),to=max(log10(rowMeans(dge$counts))),by=((max(log10(rowMeans(dge$counts)))-min(log10(rowMeans(dge$counts))))/(a-1))))
dev.off()
##----
##Cumulative percent of average counts
##----

##split different groups
groups <- levels(as.factor(metadata$time))
which(as.factor(metadata$time)==groups[1])
which(as.factor(metadata$time)==groups[2])

##rarefaction
rareify <- t(dge$counts)
rareify <- round(rareify)
raremax <- min(rowSums(rareify))
jpeg(filename = "diagnosticPlots/rarefaction.jpeg")
rarecurve(rareify, step = 20, sample = raremax, col = "blue", cex = 0.6)
dev.off()
rarefy(rareify, sample = raremax)
# plot(sum)

###------
##rarefaction
###------

##prim counter
k <- 0
percentseq <- seq(2.5,100,2.5)
cumcontr <- matrix(0,ncol = ncol(dge$counts),nrow = length(percentseq))
##loop through percent steps
for(i in percentseq) 
    { k <- k+1
##loop through samples
  for(j in 1:ncol(dge$counts))
      {
      percent <- sort(dge$counts[,j],decreasing = T)
      maxvalue <- max(percent)
      totalvalue <- sum(percent)
      length(percent)
      colnum <- (i/100)*length(percent)
      totalpercent <- (sum(percent[1:colnum])/totalvalue)*100
      cumcontr[k,j] <- totalpercent
      }
  
    }

jpeg(filename = "diagnosticPlots/cumulative_plot.jpeg")
plot(percentseq,cumcontr[,1],xlab="% of genes contributing",ylab="% of total count",type="l",col="red",main="cumulative plot")
for(j in 2:ncol(dge$counts))
  {
lines(percentseq,cumcontr[,j],col="black",lty=j)
  
}
dev.off()



## Define contrasts. ************** MODIFY ************** 
(contrasts <- as.data.frame(makeContrasts(a2vsa5=timea5, levels = des)))

## Perform tests
signif3 <- function(x) signif(x, digits = 3)
edgeR_res <- lapply(contrasts, function(cm) {
  qlf <- glmQLFTest(qlfit, contrast = cm)
  tt <- topTags(qlf, n = Inf, sort.by = "none")$table
  tt %>% dplyr::mutate_if(is.numeric, signif3)
})


## Write results to text files
if (class(edgeR_res) == "data.frame") {
  write.table(edgeR_res %>% dplyr::arrange(PValue), 
              file = gsub("rds$", "txt", outrds), 
              sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
} else {
  for (nm in names(edgeR_res)) {
    write.table(edgeR_res[[nm]] %>% dplyr::arrange(PValue), 
                file = gsub("\\.rds$", paste0("_", nm, ".txt"), outrds), 
                sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  }
}
saveRDS(list(results = edgeR_res, data = dge0), file = outrds)


sessionInfo()
date()

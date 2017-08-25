source("http://www.bioconductor.org/biocLite.R")
#biocLite('affyPLM')
#biocLite('limma')
#biocLite('GEOquery')
#biocLite('vsn')
#biocLite("makecdfenv")
library(makecdfenv)
library('affyPLM')
library('limma')
library('GEOquery')
library('R.utils')
library('vsn')

targets<-readTargets('Targets.txt')
filePaths<-paste(targets$FileNam); #wybiera kolumne "FileName" z pliku targets

for (i in 1:length(filePaths)){
	get=getGEOSuppFiles(filePaths[i],makeDirectory = FALSE, baseDir =('C:/downloads'))
}
files=(Sys.glob('*.gz'))
for (i in 1:length(files)){
  unpack=gunzip(files[[i]], overwrite=TRUE,remove=TRUE)
  print(paste0('Wypakowano',' ', i, ' ', 'elementów na ', length(files)))
}

abatch=ReadAffy();
abatch.cdf=cdfName(abatch); cat("CDF:",abatch.cdf,"\n")

abatch.gc=bg.adjust.gcrma(abatch, affinity.source="reference", type="affinities",
                           GSB.adjust=T, fast=F, optical.correct=F);

#Normalizacja
set.seed(20130517);  # dla powtarzalnosci próbkowania
abatch.vsn=normalize.AffyBatch.vsn(abatch.gc,lts.quantile=0.6,subsample=1000,
                              	log2asymp=T,log2scale=F);

pset= fitPLM(abatch.vsn,background=F,normalize=F)
e=coefs(pset);
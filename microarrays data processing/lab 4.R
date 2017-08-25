source("http://www.bioconductor.org/biocLite.R");

#biocLite('affyPLM');
#biocLite('MASS');
#biocLite('limma');
#biocLite('topGO');
#biocLite('Rgraphviz');
#biocLite('hgu133plus2.db');
#biocLite('makecdfenv');
#biocLite("GEOquery");
#biocLite("R.utils");

library(MASS)
library(vsn)
library(limma)
library(topGO)
library(Rgraphviz)
library(hgu133plus2.db)
library(makecdfenv)
library(affyPLM)
library(GEOquery)
library(R.utils)
library(vsn)

targets<-readTargets('Targets.txt')

filePaths<-paste(targets$FileNam)

for (i in 1:length(filePaths)){
	get=getGEOSuppFiles(filePaths[i],makeDirectory = FALSE, basedir='c:/downloads'))
	print(paste0('Pobrano ', i ,' na ', length(filePaths), 'elementów'))
}


files=(Sys.glob('*.gz'))
for (i in 1:length(files))
{
  unpack=gunzip(files[[i]], overwrite=TRUE,remove=TRUE)
  print(paste0('Wypakowano',' ', i, ' ', 'elementów na ', length(files)))
}

abatch=ReadAffy()

abatch.cdf=cdfName(abatch); cat("CDF:",abatch.cdf,"\n")

abatch.gc=bg.adjust.gcrma(abatch, affinity.source="reference", type="affinities", GSB.adjust=T, fast=F, optical.correct=F)

abatch.raw=abatch


set.seed(20130517)

abatch.vsn=normalize.AffyBatch.vsn(abatch.gc,lts.quantile=0.6,subsample=1000, log2asymp=T,log2scale=F)

pset= fitPLM(abatch.vsn,background=F,normalize=F)

e=coefs(pset)

for(stress in unique(targets$stress))
{
    X11();
    par(mfrow=c(2,2))
    pairs(e[,targets$stress==stress], pch=".", main=stress)
}
for(rep in unique(targets$Rep))
{
    X11();
    par(mfrow=c(2,2))
    pairs(e[,targets$Rep==rep], pch=".", main=paste("replicate", rep))
}


levels=sort(unique(targets$stress))

L.=factor(targets$stress,levels=c('medium',levels[levels!='medium']))

design=model.matrix(~ -1 + L.)
rownames(design)=L.

contrasts.matrix=makeContrasts(base  = L.medium,
                                virus = L.RV16 - L.medium,
                                smoke = L.CSE  - L.medium,
                                xaction= (  ( L.RV16plusCSE - L.medium )
                                - ( L.RV16 - L.medium + L.CSE - L.medium ) ),
                                levels = design)

								
fit=lmFit(e,design)

fit.eb=eBayes(fit)

fit2=contrasts.fit(fit,contrasts.matrix);
fit2.eb=eBayes(fit2)


adj.method="holm";
p.val = 0.05 
cont = 'virus'

tTable=topTable(fit2.eb,adjust=adj.method,coef=cont,sort.by="P",number=Inf)
sum.max=sum(tTable$adj.P.Val < p.val)
cat(cont," ", sum.max,"\n")
if(sum.max >0)
{
    print(tTable[1:sum.max,])
}

topDiffThresh=0.05;
topDiffGenes = function (score,threshold=topDiffThresh)
{
  return(score < threshold)
}

allGenes = tTable$adj.P.Val
names(allGenes) = rownames(tTable)
affyLib = paste(annotation(abatch), "db", sep = ".")

sampleGOdata = new("topGOdata", 
                    description = "Simple session", ontology = "BP",
                    allGenes = allGenes, geneSel = topDiffGenes,
                    nodeSize = 10,
                    annot = annFUN.db, affyLib = affyLib)

resultFisher = runTest(sampleGOdata, algorithm = "classic", statistic = "fisher")

rTable = GenTable(sampleGOdata, classicFisher = resultFisher,
                   orderBy = "classicFisher", topNodes = 20)
cat("sigGenes:",length(sigGenes(sampleGOdata)),"\n")



pdf("przyklad.pdf")

showSigOfNodes(sampleGOdata, score(resultFisher), firstSigNodes = 5, useInfo = 'all')

dev.off()

print(paste0('Koniec skryptu'))
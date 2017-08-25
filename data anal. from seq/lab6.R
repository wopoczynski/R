source("https://bioconductor.org/biocLite.R");
#biocLite("seqLogo");
library(seqLogo);
library(ShortRead);
library(seqinr);

r=readFastq("ERR237803.fastq.gz")
s=sread(r)
M=matrix(nrow=4, ncol=36)
M=alphabetByCycle(s,as.prob=TRUE)[1:4,]/colSums(alphabetByCycle(s,as.prob=TRUE)[1:4,])
M=makePWM(M);
seqLogo(M,ic.scale=FALSE)
test=M@consensus;

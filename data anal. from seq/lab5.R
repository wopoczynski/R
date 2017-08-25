source("http://bioconductor.org/biocLite.R");
#biocLite("ShortRead");
#biocLite("seqinr");
library(ShortRead);
library(seqinr);
seq=readFastq('ERR237803.fastq.gz');
read=sread(seq);
qual=quality(seq);
N= 338264;L= 36;G=4.6*1000000;
a=(N*L)/G;
str=toString(qual@quality);
surowy=charToRaw(str);
val=strtoi(surowy, base=16L);
finalval=val-32;
encode=encoding(qual);
###################################
probk=sample(finalval,100000,replace=FALSE);
boxplot(probk,outline=FALSE,col='green');
plot(probk))
procent=alphabetByCycle(qual)
################################################

source("https://bioconductor.org/biocLite.R");
#biocLite("seqLogo");
library(seqLogo);
library(ShortRead);
library(seqinr);
seq=readFastq('ERR237803.fastq.gz');
read=sread(seq);
dane=matrix(nrow=4, ncol=36)
dane=alphabetByCycle(read,c('A','C','G','T'));
pwm=makePWM(dane, as.prob=T);
seqLogo(pwm, ic.scale=FALSE, xaxis=TRUE, yaxis=TRUE, xfontsize=15, yfontsize=15)
r=readFastq("ERR237803.fastq.gz")
s=sread(r)
M=matrix(nrow=4, ncol=36)
M=alphabetByCycle(s,as.prob=TRUE)[1:4,]/colSums(alphabetByCycle(s,as.prob=TRUE)[1:4,])
M=makePWM(M);
seqLogo(M,ic.scale=FALSE)
wykres=c('G','T','T','C','A','G','A','A','C','T','G','C','G','C','G','G','A','A','G','G','G','G','G','G','G','G','G','G','G','A','G','G','G','G','G','G');
wykrespwm=makePWM(wykres);
seqLogo(primerpwm);
test=M@consensus;
len=length(N)
A=alphabetByCycle(s)['A',]
G=alphabetByCycle(s)['G',]
T=alphabetByCycle(s)['T',]
C=alphabetByCycle(s)['C',]
library(seqLogo);library(ShortRead);library(seqinr);
r=readFastq("ERR237803.fastq.gz");
s=sread(r);
M=matrix(nrow=4, ncol=36);
M=alphabetByCycle(s,as.prob=TRUE)[1:4,]/colSums(alphabetByCycle(s,as.prob=TRUE)[1:4,]);
M=makePWM(M);
seqLogo(M,ic.scale=FALSE);
test=M@consensus;
N=alphabetByCycle(s)['N',];
len=length(N);
A=alphabetByCycle(s)['A',];
G=alphabetByCycle(s)['G',];
T=alphabetByCycle(s)['T',];
C=alphabetByCycle(s)['C',];

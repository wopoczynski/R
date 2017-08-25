source('http://bioconductor.org/biocLite.R');
#biocLite('seqinr');
#biocLite('Biostrings');
library('seqinr');
library('Biostrings');
szympans=readDNAStringSet('szympans.fasta');
mysz=readDNAStringSet('mysz domowa.fasta');
czlowiek=readDNAStringSet('czlowiek.fasta');
p1=c();
p2=c();
p3=c();
for (i in 1:3){
  p1=c(p1, pairwiseAlignment(czlowiek[i],mysz[i], type='global', gapOpening=-5, scoreOnly = T));
  p2=c(p2, pairwiseAlignment(czlowiek[i],szympans[i], type='global', gapOpening=-5, scoreOnly = T));
  p3=c(p3, pairwiseAlignment(mysz[i],szympans[i], type='global', gapOpening=-5, scoreOnly = T));
}
dopasowanie=matrix(nrow=3,ncol=3);
gen1=as.matrix(c(p1[1], p2[1], p3[1]));
gen2=as.matrix(c(p1[2], p2[2], p3[2]));
gen3=as.matrix(c(p1[3], p2[3], p3[3]));
dopasowanie=cbind(gen1,gen2,gen3);
colnames(dopasowanie)=c('MB','MC1R','HBB');
rownames(dopasowanie) = c('Czlowiek - Mysz','Czlowiek - Szympans','Mysz - Szympans')
czlowiekproc=alphabetFrequency(czlowiek,as.prob = T);
myszproc=alphabetFrequency(mysz,as.prob = T);
szympansproc=alphabetFrequency(szympans,as.prob = T);

nukl=c('A','T','G','C');
procgen1=cbind(czlowiekproc[1,nukl],myszproc[1,nukl],szympansproc[1,nukl]);
colnames(procgen1)=c('czlowiek','mysz','szympans');
procgen2=cbind(czlowiekproc[2,nukl],myszproc[2,nukl],szympansproc[2,nukl]);
colnames(procgen2)=c('czlowiek','mysz','szympans');
procgen3=cbind(czlowiekproc[3,nukl],myszproc[3,nukl],szympansproc[3,nukl]);
colnames(procgen3)=c('czlowiek','mysz','szympans');
png(filename="1gen.png")
barplot(procgen1, main="Udzia³ procentowy Nukleotydów w MB ", ylab= "Procent",beside=T, col=c("black","darkblue","darkgreen","darkgrey"));
legend("topleft", c("A","T","G","C"), cex=1, bty="n", fill=c("black","darkblue","darkgreen","darkgrey"));
dev.off()

png(filename="2gen.png")
barplot(procgen2, main="Udzia³ procentowy Nukleotydów w MC1R ", ylab= "Procent",beside=T, col=c("black","darkblue","darkgreen","darkgrey"));
legend("topleft", c("A","T","G","C"), cex=1, bty="n", fill=c("black","darkblue","darkgreen","darkgrey"));
dev.off()

png(filename="3gen.png")
barplot(procgen3, main="Udzia³ procentowy Nukleotydów w HBB ", ylab= "Procent",beside=T, col=c("black","darkblue","darkgreen","darkgrey"));
legend("topleft", c("A","T","G","C"), cex=1, bty="n", fill=c("black","darkblue","darkgreen","darkgrey"));
dev.off()

szympans1=read.fasta('szympans.fasta');
mysz1=read.fasta('mysz domowa.fasta');
human1=read.fasta('czlowiek.fasta');
humanaat=count(human1[[1]],3);
humanaat2=humanaat['aat'];
myszaat=count(mysz1[[1]],3);
myszaat2=myszaat['aat'];
szympansaat=count(szympans1[[1]],3);
szympansaat2=szympansaat['aat'];
tabaat1gen=cbind(humanaat2,myszaat2,szympansaat2);
colnames(tabaat1gen)=c('H. sapiens','M. musculus', 'P. troglodytes');

humanaat=count(human1[[2]],3);
humanaat2=humanaat['aat'];
myszaat=count(mysz1[[2]],3);
myszaat2=myszaat['aat'];
szympansaat=count(szympans1[[2]],3);
szympansaat2=szympansaat['aat'];
tabaat2gen=cbind(humanaat2,myszaat2,szympansaat2);
colnames(tabaat2gen)=c('H. sapiens','M. musculus', 'P. troglodytes');

humanaat=count(human1[[3]],3);
humanaat2=humanaat['aat'];
myszaat=count(mysz1[[3]],3);
myszaat2=myszaat['aat'];
szympansaat=count(szympans1[[3]],3);
szympansaat2=szympansaat['aat'];
tabaat3gen=cbind(humanaat2,myszaat2,szympansaat2);
colnames(tabaat3gen)=c('H. sapiens','M. musculus', 'P. troglodytes');


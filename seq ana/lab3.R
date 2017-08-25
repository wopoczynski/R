ptandem=read.table('wynik.txt', header = T)


report=read.table("wynik.txt", header = T)
length=report[['Length']]
which.max(length)
report[193,]

ratio=report[[4]]
which.max(ratio)
report[37,]


hist(length, breaks=(min(length)-0.5):(max(length)+0.5),)
sum(report[,3]==2)
table(report[,3])

dwojki=report[which(report[,3]==2),]
trojki=report[which(report[,3]==3),]

dwojki_w=unique(dwojki[c("Motif")])
trojki_w=unique(trojki[c("Motif")])




write.table(dwojki_w,file="dwojki.txt")
write.table(trojki_w,file="trojki.txt")
r=read.table("dwojki.txt")
w=read.table("trojki.txt")

l4=c();
l3=c();
l2=c();



report_1=read.table("seq1.txt", header = T)
l4= c(l4,sum(report_1[,3]==4))
l3= c(l3,sum(report_1[,3]==3))
l2= c(l2,sum(report_1[,3]==2))
report_2=read.table("seq2.txt", header = T)
l4= c(l4,sum(report_2[,3]==4))
l3= c(l3,sum(report_2[,3]==3))
l2= c(l2,sum(report_2[,3]==2))
report_3=read.table("seq3.txt", header = T)
l4= c(l4,sum(report_3[,3]==4))
l3= c(l3,sum(report_3[,3]==3))
l2= c(l2,sum(report_3[,3]==2))
report_4=read.table("seq4.txt", header = T)
l4= c(l4,sum(report_4[,3]==4))
l3= c(l3,sum(report_4[,3]==3))
l2= c(l2,sum(report_4[,3]==2))
report_5=read.table("seq5.txt", header = T)
l4= c(l4,sum(report_5[,3]==4))
l3= c(l3,sum(report_5[,3]==3))
l2= c(l2,sum(report_5[,3]==2))

genle=c(2080, 9618, 39852, 103343, 195834);

plot(genle,l2,type='b',xlab='Dlugoœæ genomu',ylab='Iloœæ sekwencji tandemowych',col='red',ylim=c(0.0,7000.0))
par(new=T)
plot(genle,l3,type='b',axes=F,xlab='',ylab='',col='blue',ylim=c(0.0,7000.0))
par(new=T)
plot(genle,l4,type='b',axes=F,xlab='',ylab='',col='black',ylim=c(0.0,7000.0))
title('Zaleznosc ilosci sekwencji tandemowych od dlugosci genomu')
legend('topleft', c('dlugosc: 2','dlugosc:3','dlugosc:4'), lty=c(1,1), col=c('red', 'blue', 'black'), bty='n', cex=.75)




lsp=c();
l3p=c();
l2p=c();
report_1p=read.table("seq1p.txt", header = T)
lsp= c(lsp,sum(report_1p[,3]))
l3p= c(l3p,sum(report_1p[,3]==3))
l2p= c(l2p,sum(report_1p[,3]==2))
report_2p=read.table("seq2p.txt", header = T)
lsp= c(lsp,sum(report_2p[,3]))
l3p= c(l3p,sum(report_2p[,3]==3))
l2p= c(l2p,sum(report_2p[,3]==2))
report_3p=read.table("seq3p.txt", header = T)
lsp= c(lsp,sum(report_3p[,3]))
l3p= c(l3p,sum(report_3p[,3]==3))
l2p= c(l2p,sum(report_3p[,3]==2))
report_4p=read.table("seq4p.txt", header = T)
lsp= c(lsp,sum(report_4p[,3]))
l3p= c(l3p,sum(report_4p[,3]==3))
l2p= c(l2p,sum(report_4p[,3]==2))
report_5p=read.table("seq5p.txt", header = T)
lsp= c(lsp,sum(report_5p[,3]))
l3p= c(l3p,sum(report_5p[,3]==3))
l2p= c(l2p,sum(report_5p[,3]==2))
report_6p=read.table("seq6p.txt", header = T)
lsp= c(lsp,sum(report_5p[,3]))
l3p= c(l3p,sum(report_5p[,3]==3))
l2p= c(l2p,sum(report_5p[,3]==2))

LAt=matrix(lsp, nrow=6, ncol=1)
L3t=matrix(l3p, nrow=6, ncol=1)
L2t=matrix(l2p, nrow=6, ncol=1)

rownames(LAt) <- c('H. Sapiens', 'P. Troglodytes', 'C. Lupus', 'M. Musculus', 'D. Rerio','D. Melanogaster')
colnames(LAt) <- c('Wszystkie')
colnames(L3t) <- c('3')
colnames(L2t) <- c('2')
LT=cbind(LAt,L3t,L2t)
LT2=t(LT)
barplot(LT2, main = "Iloœæ sekwencji tandemowych",  xlab = "Organizm",  col = c("blue","red","black"),  beside=TRUE)
legend("topright", c("Wszystkie","dlugosc: 3","dlugosc: 2"), fill = c("blue","red","black"))
A=c(1,2,3,4,5,6)
B=matrix(A, byrow=T,nrow=3) 
print(B) 
summary(B)
suma=sum(B) 
C=matrix(seq(6),byrow=F,nrow=3) 
print(C) 
D=cbind(B,C) 
colnames(D)=c('Darek','Darek','Anna','Jarek') 
rownames(D)=c('waga','wiek','wzrost') 
sum(D['waga',])
D=D[,-2] 
czas=c(20,25,10) 
D=rbind(D,czas) 
dane=data.frame(t(D)) 
stan=matrix(c('chory','chory','zdrowy')) 
dane=data.frame(dane,stan) 
dane['Darek','stan'] 
chorzy=which(dane[,'stan']=='chory') 
index=sort(dane$wiek,index.return=T) 
setwd('C:/Users/wojte/Dropbox/studia/GFiP/lab1');
dane1=dane[chorzy,] 
dane1['waga'>2] 
write.table(file='dane.txt',dane1, sep=';') 
dane2=read.table('dane.txt',sep=';')


data(ChickWeight)
Weight=ChickWeight[,'weight']
Time=ChickWeight[,'Time']
Chick=ChickWeight[,'Chick']
Diet=ChickWeight[,'Diet']
index=sort(ChickWeight$Time,index.return = T)
dane=ChickWeight[index$ix,]
result=head(dane,n=20)
Weight1=result[,'weight']
Time1=result[,'Time']
Chick1=result[,'Chick']
Diet1=result[,'Diet']
boxplot(Weight1,col='red',xlab="wartosci",main="Wykres", outline = F)#boxplot bez punktóW wykraczaj¹cych


n=length(Weight1);
sdW1=sd(Weight1)
sredniawaga1=mean(Weight1)
#przedzial ufnosci 95%
przedzialufnosci=round(sredniawaga1+c(-1,1)*sdW1/sqrt(n)*qnorm(.975),2)


t.test(Weight, Weight1)

write.table(file='dane.txt',result, sep=';')
library("seqinr")
d=read.fasta("d.fasta")
m=read.fasta("m.fasta")
d_seq=d[[1]]
m_seq=m[[1]]
d_dl=length(d_seq)
m_dl=length(m_seq)

a=1
ramka=100
b=ramka
zawartoscGC=c()
di=floor(d_dl/ramka)
skok_okna=c(1:di)
for (i in 1:di){
	frag=d_seq[a:ramka]
	y=GC(frag)
	zawartoscGC[i]=y
	a=a+ramka
	b=b+ramka
}
plot(skok_okna,zawartoscGC, "l")

myc = read.fasta("m.fasta")
myc_seq = myc[[1]]
myc_dl = length(myc_seq)
ramka = 20000
a = 1
b = ramka
y = c()
skok = floor(myc_dl/ramka)
skok_vec = c(1:skok)
for (i in 1:skok){
	fragment = myc_seq[a:b]
	zaw_gc = GC(fragment)
	y[i] = zaw_gc
	a = a+ramka
	b = b+ramka	
}
plot(skok_vec, y, "l")

zawAT = function(sekw){
	zasady=table(sekw)
	at=zasady[["a"]]+zasady[["t"]]
}
zliczAT=zawAT(m)
print(zliczAT)

ramka=20000
a=1
b=ramka
y=c()
di=floor(m_dl/ramka)
skok_okna=c(1:di)
for (i in 1:di){
	frag=m_seq[a:b]
	zasady=table(frag)
	iloscat=zasady[['a']]+zasady[['t']]
	udzial=iloscat/ramka
	y[i]=udzial
	a=a+ramka
	b=b+ramka
}
plot(skok_okna,y, "l")

library(seqinr)

myc=read.fasta("m.fasta")
myc_seq=myc[[1]]

myc_len=length(myc_seq)

bases=table(myc_seq)
perc_g=bases[['g']]/myc_len
perc_a=bases[['a']]/myc_len
perc_c=bases[['c']]/myc_len

prob_gac=perc_g*perc_a*perc_c

trimers=count(myc_seq, 3)
quant_gac=trimers[["gac"]]
perc_gac=quant_gac/(myc_len-2)

if (perc_gac>prob_gac){
	print("S³owo 'GAC' jest nadreprezentowane")
} else if (perc_gac < prob_gac){
	print("S³owo 'GAC' jest niedoreprezentowane")
} else{
	print("S³owo 'GAC' odpowiada wyliczonemu prawdopodobieñstwu")
}
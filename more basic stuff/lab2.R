#source('http://bioconductor.org/biocLite.R');
#biocLite('Biostrings');
#biocLite('hgu95av2.db');
#biocLite('annotate');
#library('Biostrings');
#library('hgu95av2.db');
#library('annotate');

#hgu95av2ACCNUM
plik=read.table('affy.txt');
tab=hgu95av2ACCNUM;
mapped_probes=mappedkeys(tab);
tab_lista=as.list(tab[mapped_probes]);
if(length(tab_lista) > 0) {
  tab_lista_w=tab_lista[1:5];
}

#getSEQ
wektorseq=c();
for (i in 1:length(tab_lista_w)){
  wektorseq[i]=getSEQ(tab_lista[[i]])
}


#pairwiseAlignmen
wynik=c()
wyn=matrix(ncol=5, nrow=5)
for (i in 1:5){
  for (j in 1:5){
    if (i>j){
      wynik=pairwiseAlignment(wektorseq[i],wektorseq[j], type='global', gapOpening=-5)
      wyn[i,j]=wynik@score
    }
  }
}
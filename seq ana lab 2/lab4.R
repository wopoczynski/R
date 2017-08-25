dlchrom1=234560296;
dlchrom3=203221608;
dlchrom4=173516952;
dlchrom8=150392723;
dlchrom12=108540874;
dlchrom18=75010233;
contig1=1+2698+161+14809;
contig3=1+2323+171+11569;
contig4=1+9831+71+1664;
contig8=1215+1+7950+63;
contig12=820+1+28+5707;
contig18=690+1+3904+39

dlugosc=c(dlchrom1,dlchrom3,dlchrom4,dlchrom8,dlchrom12,dlchrom18);
contig=c(contig1,contig3,contig4,contig8,contig12,contig18);

plot(dlugosc,contig, type="b", xlab="dlugosc sekwencji", ylab="ilosc contigow",);
title("zaleznosci ilosci contigow od dlugosci sekwencji chromosomu");


lspokrycie=48.0;
lscontig=458935;
penpokrycie=61;
pencontig=230930;
varpokrycie=81;
varcontig=110959;
halpokrycie=103;
halcontig=31786;
plasmidpokrycie=136;
plasmidcontig=349;
simpokrycie=180;
simcontig=82;
pokrycie=c(lspokrycie,penpokrycie,varpokrycie,halpokrycie,plasmidpokrycie,simpokrycie);
contigC=c(lscontig,pencontig,varcontig,halcontig,plasmidcontig,simcontig);
plot(pokrycie,contigC,type='b',xlab='pokrycie genomu',ylab='ilosc contigow');
title("zaleznosc ilosci contigow od pokrycia genomu");

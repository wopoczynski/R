library(Biostrings)
library(seqinr)
library(BiocGenerics)

primer_prawy = ("AAGCCAACTTCAACACCCCA");
primer_lewy = ("GAAGGTGGAGGCTGACATCC");

primer_prawy_split = strsplit(primer_prawy ,split="");
primer_lewy_split = strsplit(primer_lewy,split="");

ile_nukleotydow_p = table(primer_prawy_split);
ile_nukleotydow_l = table(primer_lewy_split);

temperatura_topnienia_p = 2*(ile_nukleotydow_p["A"]+ile_nukleotydow_p["T"])+4*(ile_nukleotydow_p["G"]+ile_nukleotydow_p["C"]);#60st
temperatura_topnienia_l = 2*(ile_nukleotydow_l["A"]+ile_nukleotydow_l["T"])+4*(ile_nukleotydow_l["G"]+ile_nukleotydow_l["C"]);#64st

zawartoscGC_p = (ile_nukleotydow_p["G"]+ile_nukleotydow_p["C"])/(sum(ile_nukleotydow_p))*100; #50%
zawartoscGC_l = (ile_nukleotydow_l["G"]+ile_nukleotydow_l["C"])/(sum(ile_nukleotydow_l))*100; #60%

ostatnie_primer_p = head(primer_prawy_split[[1]], n=5);#3' POCZATKOWY ## AAGCC
ostatnie_primer_l = tail(primer_lewy_split[[1]], n=5);#3' KONCOWY ## CATCC

odwrotny_koniec_primer_p = rev(ostatnie_primer_p[]);
odwrotny_koniec_primer_l =rev(ostatnie_primer_l[]);
komplementarnosc_p = pairwiseAlignment(ostatnie_primer_p[],odwrotny_koniec_primer_p[],type="global");
komplementarnosc_l = pairwiseAlignment(ostatnie_primer_l[],odwrotny_koniec_primer_p[],type="global");

wzorce_poliX = c("AAAAA","GGGGG","TTTTT","CCCCC");
poliX_p = c();
poliX_l = c();

for (i in 1:length(primer_prawy[[1]]))
{
  poliX_p = grepl(wzorce_poliX[i],primer_prawy);
  poliX_l = grepl(wzorce_poliX[i],primer_lewy);
}

#powtarzajace sie sekwencje
wzorce_powtorzenia = c("ATATATATAT","AGAGAGAGAG","ACACACACAC",
         "TATAATATAT","TCTCTCTCTC","TGTGTGTGTG",
         "CTCTCTCTCT","CGCGCGCGCG","CACACCACAC",
         "GCGCGCGCGC","GAGAGAGAGA","GTGTGTGTGT");
powtorzenia_p = c();
powtorzenia_l = c();

for (i in 1:length(primer_prawy[[1]]))
     {
  powtorzenia_p = grepl(wzorce_powtorzenia[i],primer_prawy[1]);
  powtorzenia_l = grepl(wzorce_powtorzenia[i],primer_lewy[1]);
}

source("http://bioconductor.org/biocLite.R")
biocLite("ShortRead")
biocLite("seqinr")
biocLite("seqLogo")
library(ShortRead)
library(seqinr)
library(seqLogo)



seq = readFastq('ERR237803.fastq.gz')

read = sread(seq)
qual = quality(seq)

#a = NL/G
#L - dlugosc odczytu
#N - ilosc odczytow
#G - dlugosc genomu


a = 338264*36/4600000


str = toString(qual@quality)

str2 = charToRaw(str)

val = strtoi(str2, base = 16L)
val2 = val-32

encode = encoding(qual)


samp_seq = sample(seq, 100000, replace = F)
samp_read = sread(samp_seq)
samp_qual = quality(samp_seq)

samp_vals = matrix(nrow = 100000, ncol = 36)
for(i in 1:100000){
  vals_str = toString(samp_qual@quality[i])
  vals_str2 = charToRaw(vals_str)
  samp_vals[i,] = strtoi(vals_str2, base = 16L)-32
}


boxplot(samp_vals, outline =  FALSE, col = 'green',
        main = 'Wykres pudelkowy wskaznikow jakosci dla danej pozycji odczytu',
        xlab = 'Pozycja odczytu', ylab = 'Wartosc wskaznika jakosci')

		lines(colMeans(samp_vals, na.rm = TRUE), col = 'red', type = 'l')

matrix_cyc = matrix(nrow = 4, ncol = 36)

matrix_cyc = alphabetByCycle(samp_read)[1:4,]/colSums(alphabetByCycle(samp_read))

matrix_cyc_per = matrix_cyc*100

matrix_cyc_per_tr = t(matrix_cyc_per)


matplot(c(1:36) ,matrix_cyc_per_tr, pch = 1, type = 'o', col=c("red","green","blue", "yellow"),
        main = 'Procentowa zawartosc nukleotydow dla danej pozycji odczytu',
        xlab = 'Pozycja odczytu', ylab = 'Procentowa zawartość')
legend("topright", legend = c('A', 'C', 'G', 'T'),
       col = c("red", "green", "blue", "yellow"), cex = 0.8, lty = 1)


	   matrix_cyc_N = alphabetByCycle(samp_read)[15,]/colSums(alphabetByCycle(samp_read))

	   matrix_cyc_N_per = matrix_cyc_N*100


	   matplot(c(1:36), matrix_cyc_N_per, pch = 1, type = 'o',
        main = 'Procentowa zawartosc nierozpoznanych nukleotydow dla danej pozycji odczytu',
        xlab = 'Pozycja odczytu', ylab = 'Procentowa zawartość')


		CG_cyc_per = letterFrequency(samp_read, letters = 'GC', as.prob = TRUE)

		CG_cyc_per = CG_cyc_per*100


		hist(CG_cyc_per, breaks = 40, main = 'Procentowa zawartosc par GC na danej pozycji odczytu',
     xlab = 'Procentowa zawartosc', ylab = 'Ilosc zliczen')


	 

	 pwm_cyc = makePWM(matrix_cyc)

	 seqLogo(pwm_cyc, ic.scale = FALSE)


	 matrix_cyc_count = alphabetByCycle(samp_read)[1:4,]


	 logo_high = c('G', 'T', 'T', 'C', 'A', 'G', 'A', 'A', 
         'C', 'T', 'G', 'C', 'G', 'C', 'G', 'G', 'A', 'A',
         'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G',
         'A', 'G', 'G', 'G', 'G', 'G', 'G')


		 logo_high_num = c(3, 4, 4, 2, 1, 3, 1, 1, 2, 4, 3, 2, 3, 2,
                  3, 3, 1, 1, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
                  1, 3, 3, 3, 3, 3, 3)

				  test_count = pwm_cyc@consensus


				  

				  N_vec = alphabetByCycle(samp_read)[15, ]

				  for (i in 1:36){
  matrix_cyc_count[logo_high_num[i], i] = matrix_cyc_count[logo_high_num[i], i] + N_vec[i]
}

korekta = matrix_cyc_count/colSums(matrix_cyc_count)
pwm_korekta = makePWM(korekta)

seqLogo(pwm_korekta, ic.scale = FALSE)

korekta_per = korekta*100
korekta_per_tr = t(korekta_per)


matplot(c(1:36) ,korekta_per_tr, pch = 1, type = 'o', col=c("red","green","blue", "yellow"),
        main = 'Procentowa zawartosc nukleotydow dla danej pozycji odczytu po korekcie',
        xlab = 'Pozycja odczytu', ylab = 'Procentowa zawartosc')
legend("topright", legend = c('A', 'C', 'G', 'T'),
       col = c("red", "green", "blue", "yellow"), cex = 0.8, lty = 1)


	   samp_mat = as.matrix(samp_read)

	   for(i in 1:100000){
  for(j in 1:36){
    if(samp_mat[i, j] == "N"){
      samp_mat[i, j] = logo_high[j]
    }
  }
}  
#

CG_cyc_per_kor = letterFrequency(samp_mat, letters = 'GC', as.prob = TRUE)

CG_cyc_per_kor = CG_cyc_per_kor*100


hist(CG_cyc_per_kor, breaks = 40, main = 'Procentowa zawartosc par GC na danej pozycji odczytu po korekcie',
     xlab = 'Procentowa zawartosc', ylab = 'Ilosc zliczen')

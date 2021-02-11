#  R codes
# Author: Hui-Heng Lin, PhD. Copyright reserved. 



""" On 6th may 2020 @ rstudio of winserver 2016 at yuebei hopistal's 3-13f-1341 office room; created since 6th may 2020, continuous updates, modifications presented

Previously stupidly forgotten to save all data and variables in R, now i have to reconstruct the svm model from zero beginning!

above reconstruction is Fine, but this time, has to make it a clean and standard code!!!!======hehe still failed, but make it on 11th may 2020@ notepad ++ win8.1 on HP omen 17

## remember to set.seed(x) stupid for reproducibility!

Author: Hui-heng Lin , phd. postdoc   



R version 3.6.3 (2020-02-29) -- "Holding the Windsock"
Copyright (C) 2020 The R Foundation for Statistical Computing
Platform: x86_64-w64-mingw32/x64 (64-bit)

R是自由软件，不带任何担保。在某些条件下你可以将其自由散布。用'license()'或'licence()'来看散布的详细条件。

R是个合作计划，有许多人为之做出了贡献. 用'contributors()'来看合作者的详细情况. 用'citation()'会告诉你如何在出版物中正确地引用R或R程序包。

用'demo()'来看一些示范程序，用'help()'来阅读在线帮助文件，或 用'help.start()'通过HTML浏览器来看帮助文件。 用'q()'退出R.


> citation()

To cite R in publications use:   R Core Team (2020). R: A language and  environment for statistical computing. R   Foundation for Statistical Computing, Vienna,   Austria. URL https://www.R-project.org/. ========odel R-project's base was in Vienna?

A BibTeX entry for LaTeX users is

  @Manual{,
    title = {R: A Language and Environment for Statistical Computing},
    author = {{R Core Team}},
    organization = {R Foundation for Statistical Computing},
    address = {Vienna, Austria},
    year = {2020},
    url = {https://www.R-project.org/},
  }

We have invested a lot of time and effort in
creating R, please cite it when using it for
data analysis. See also ‘citation("pkgname")’
for citing R packages.

"""

"""
this is a standard code full flow of (re-) construction of SVM, should paste index/category at the beginning of this scrip file? E.g.:

                      Index
1.0 xxxxx  line xx-xx
1.1 xxxxx  line xx-xx?
1.1.1 xxxx
1.1.2.... xxxx
1.2 xxxxx
1.3.....  xxxx


"""


# stupid R make the default working folder as c:/windows/system32!!  so stupid!!  Set in Rstudio Option or use command to change workin folder to 


setwd("D:/RworkDir/") # done. meaning set working directory



library(e1071) # load the svm lib in R
citation(e1071) # error odel!
citation("e1071") # outputed below. Correct one, require "  " on "pack_name"

"""
> citation("e1071") 

在出版物中使用程序包时引用‘e1071’:

  David Meyer, Evgenia Dimitriadou, Kurt Hornik,   Andreas Weingessel and Friedrich Leisch  (2019). e1071: Misc Functions of the
  Department of Statistics, Probability Theory   Group (Formerly: E1071), TU Wien. R package   version 1.7-3.   https://CRAN.R-project.org/package=e1071   =============hehe, lin-cih jen的libsvm要不要引用？这个interface lin 没参加开发的呵呵。

A BibTeX entry for LaTeX users is

  @Manual{,
    title = {e1071: Misc Functions of the Department of Statistics, Probability
Theory Group (Formerly: E1071), TU Wien},
    author = {David Meyer and Evgenia Dimitriadou and Kurt Hornik and Andreas Weingessel and Friedrich Leisch},
    year = {2019},
    note = {R package version 1.7-3},
    url = {https://CRAN.R-project.org/package=e1071},
  }              """

# step 1 : data preparation: sequences, vectors, labels, splitting sets: 6:2:2 as trainingset : cross-valiation set : testing set

# 1.1 get the raw/standard BRCA1 gene sequence from genome, load it from local disk file into a variable==============shit!! should make the seqeuence continuously in one-line!! multiple lines will include "\n" affecting the position calling! shit!


library(Biostrings) # need to use the DNAString() for format conversion


# try scan(). scan() was self-attached , no need to call lib nor install!
dna<-scan( file="I:/1 prjsRes/MLpredSnpVarMutants/myOwnData2020ensemblClinVarDbSPNs/fulllengthBRCA1DNAseq125969bpFromEnsemblGenomViewer.txt", what = character(0) ) 
head(dna)
class(dna) # character
dna[1] # outputed annoation:">GCh38...xxx" ==========fuck this shit, tried the shit lib 'seqinr' and 'biostring' 's shit import function, forced to import as fasta format, so I manually added " > xxxx annotation xxx" into it. In fact now, got the scan(), no need to add the shit "< xxxx annotation"

dna[2] # outputed long seq!
dna[3] # NA
dna[2][1] # shit still long seq!
dna[2][2] # shit NA! noway to get single base output []!



DNAString(dna[2]) # heyyo!!! got it!!
d3<-DNAString(dna[2]) ; class(d3) # odel "biostrings" type!
d3[2] # so cool!!!! single base letter presented eventually!
d3[1];d3[2];d3[4] # heyyo! outputed T C T

substr(d3,1,1) # T!
substr(d3,1,1) # also T! won't extract and remove!
substr(d3,2,2) # cool ! C
substr(d3,1,2) # cool ! TC!

# 1.2 copy another var for operating, do not operate the original BRCA1 seqeuence variable!

seq_op<-d3; class(seq_op) # cool, the same type: "biostrings"
seq_ts<-d3; seq_ts[5] # showed T
seq_ts[5]<-"K";
seq_ts[5] # odel! became "K"!
substr(seq_ts,5,5) # K desu
substr(seq_ts,5,5) # K still
substr(seq_ts,5,5) <-"T" # error del! 无法在非字符的对象里替换子字串 .hehe biostring 就不算string? hehe. 
# more direct way?
seq_ts[5]<-"T";
seq_ts[5] # fine be back to T now!


""" 1.2.5----- save/backup in the middle: all var into .Rdata. Backup is vital! ------"""
# can either do it via clicking graphical interface on Rstudio or use command:
save.image(file = "D:/RworkDir/reconstructSVM6may2020.RData", version = 1.0 , safe = TRUE) #  use reverse / , not the linux one :/, and the recommended suffix was ".RData", uppercase "D", not the ".Rdata". Matters?====================use load(xxx) if wanna load var/data from .RData file!




# 1.3 assign the start coordinate of BRCA1 gene in Gch38 version of human genome is 43044285
s1<-43044285; print(s1); print(s1-3) # cool


# 1.4 from  disk file load the numbers, BRCA1 mutations/SNPs location/co-ordinates, mutatated bases, label:benign/patho

df<-read.csv("I:/1 prjsRes/MLpredSnpVarMutants/myOwnData2020ensemblClinVarDbSPNs/mergedBRCA1BenignPatho1084.csv",header=T)  # I:/ indicates that data were in mobile disk!


# check 
df # cool. loaded with number, mutation location/co-ordinates, mutatated bases, label:benign/patho

str(df) # 1084 lines loaded into the data_frame type. Based on displayings, bases should be df[,4]; co-ordinates should be df[,3] 

df[580:590,2] # 1 to 586 are benign, 587 to 1086 are pathogenic! =============this was wrong!
df[586:587,2] # shit! above wrong!
df[585:587,2] # shit! it is 1 to 585 are benign. others since 586 are patho!!!!!
# so, also shit found previous trial codes were wrong(mixed patho with benign tags) in data split!
"""
bn_tag<-sample(1:583,583); bn_tag; # cool, check, no repeats ================= previous wrong shit !!  till 585 were benign!
pth_tag<-sample(584:1084,500); 
pth_tag; # cool, check, no repeats

"""


# 1.5 to mutate the seq using raw seq and bases' mutation info and location co-ordinate data. Using substr( ) function inside LOOP! No need to load library, hehe substr() was self-attached

# 1.5.1 test calling [x,y]
df[1,3] # yes, co-ordinate was it, 43094550
df[2,4] # yes, nucleotide base was it, the "C"


# 1.5.1.1 test if substr() overwrites/writeOver/change the varaible?=======tested, won't overwrite!


# 1.6  an emtpy list var is required to store each mutant of seq
mutants<-list()


# 1.7 create mutants . Two options are available: 1 save all the mutant seqs, then convert into vectors afterwards; 2 inside loop, create one mutant, vectorize right away and store into a var.=========save the middle variable, so this time choose option 1 still!

# test
df[7,3] # 43071102
seq_ts[43071102-s1] # A desu
seq_ts[df[7,3]-s1] # cool still A, hehe


# 1.7 above go!
for (i in 1:1084) {
  seq_op[df[i,3]-s1]<-df[i,4]; # replace one base each time
  mutants[i]<-seq_op; # save the mutant into a list var
  seq_op<- d3 # every time 's ending, recover it to wildtype/original standard BRCA1 seq is important!
};  dim(mutants)# ================shit!! error!!   'value' must be a DNAString object (or coercible to a DNAString object)


# fix : add as.character

for (i in 1:1084) {
  seq_op[df[i,3]-s1]<-as.character(df[i,4]); # replace one base each time
  mutants[i]<-seq_op; # save the mutant into a list var
  seq_op<- d3 # every time 's ending, recover it to wildtype/original standard BRCA1 seq is important!
};  dim(mutants) # odel NULL
mutants[1] # cool,hehe! seq shown
# to  confirm if mutations succeeded!
length(mutants) # 1084
# comparison
mutants[7]
mutants[7][3] # odel null! 
mutants[7,3] # odel error! noway to check single base? need substr()??

substr(mutants[7],43071102-s1,43071102-s1) # odel nothing !
substr(mutants[7],1,1) # odel n!
substr(mutants[7],2,3) #  odel "ew" what the fuck are these??

str(mutants[7]) # odel humor ! no seq?
mutants[7] 

mutants # list outputed
mutants[[1]] # the same with [1] a!
mutants[[1]][4] # heiyo! T outputed!!!!!!!!!!!!!!
mutants[[7]][4] # heiyo! T outputed!!!!!!!!!!!!!!
mutants[7][4] # del double embedded [[]] is required!!

substr(mutants[[7]],43071102-s1,43071102-s1) # C desu!!
d3[43071102-s1] # original base was A
df[7,4] # "C" dazo! so repalcement succeeded!! yeah!!!!




"""           1.8----- save/backup in the middle: all var into .Rdata. Backup is vital!       ------"""
# can either do it via clicking graphical interface on Rstudio or use command:

save.image(file = "D:/RworkDir/reconstructSVM6may2020.RData", version = 2.0 , safe = TRUE) #  use reverse / , not the linux one :/, and the recommended suffix was ".RData", uppercase "D", not the ".Rdata". Matters?====================use load(xxx) if wanna load var/data from .RData file!


# 1.9 convert 1084 mutant seqs into different vectors using different featurization methods, then split sample in random way!
library(BioMedR)
# 1.9.1 one-mer only 
extrDNAkmer(mutants,k=2) # error!
extrDNAkmer(as.character(mutants),k=2) # shit !! no error! but all NAs!! not recognized as ATCS
extrDNAkmer(mutants[],k=2) # error!
extrDNAkmer(as.character(mutants[]),k=2) #  # shit !! no error! but all NAs!! not recognized as ATCS
extrDNAkmer(as.character(mutants[[]]),k=2) # error!
extrDNAkmer(mutants[[1]],k=2) # error!
extrDNAkmer(as.character(mutants[[1]]),k=2) # heiyo!!!cool!!
extrDNAkmer(as.character(mutants[[1]]),k=3, upto=T) # heiyo!!!cool!!
extrDNAkmer(as.character(mutants[[2]]),k=3, upto=T) # heiyo!!!cool!!
extrDNAkmer(as.character(mutants[[3]]),k=3, upto=T) # heiyo!!!cool!!

# check  mutation points 2,3
df[2,3]; df[2,4] # 43094795, C
df[3,3]; df[3,4] # 43097246, A
d3[df[2,3]-s1] # A in orginal raw seq zo!
d3[df[3,3]-s1] # T in orginal raw seq zo!
mutants[[2]][df[2,3]-s1] # cool! C zo! mutation succeeded!
mutants[[3]][df[3,3]-s1] # cool! A zo! mutation succeeded!==========fuck shit then why  k-mer DNA all the same?

# 1.9.2 one + two + three-mers. Use "upto = T" to get all? can just extract 1-mer's feature by data_frame[x,y]?
tryvec<-extrDNAkmer(as.character(mutants[[3]]),k=3, upto=T) # heiyo!!!cool!!
str(tryvec)
class(tryvec) # integer
summary(tryvec) # heiyo outputed min. max.
length(tryvec) # 84
tryvec[1:3] # cool o! showed first 3 vectors, but with A C G!
as.numeric(tryvec[1:3]) # cool o! without A C G at head now!!


vec3mer<-list()
for (i in 1:1084){
  vec3mer[i] <- as.numeric(extrDNAkmer(as.character(mutants[[i]]), k=3, upto=T))} #  <- numeric() not work! must be "as.numeric()".  NO error but warnings!

# check
class(vec3mer) # list
str(vec3mer) #  multiple lines of $ : num 35282; $ : num 35281; $ : num 35282.  don't know what 
summary(vec3mer) # long screen
dim(vec3mer) # null
length(vec3mer) # 1084 
vec3mer[1] # odel 35281 NOT a long vector containing multiple number?
vec3mer[2] # odel 35281

vec3mer[[1]] # odel 35281!

head(vec3mer) # fuck shit!! only got the first number of the whole vector?


# again try to give double [[i]] to vec3mer storage variable!
vec3mer<-list() # clean the var again
vec3mer # fine, cleaned as empty list()
for (i in 1:10){ # firstly use 10 to test, no need to run whole!
  vec3mer[[i]] <- as.numeric(extrDNAkmer(as.character(mutants[[i]]), k=3, upto=T))} # !!!!!!!!!!!! tested confirmed !! hayali double [[i]] is required for vec3mer storage! otherwise only first number of the long vector was saved!

vec3mer # cool o! all 10 vectors were saved inside!
vec3mer[[10]] # cool, called  the 10th vector!
vec3mer[[10]][1:6] # cool, called  the first 6 numbers of 10th vector!


# fine, now run for 1084!
vec3mer<-list() # clean the var again
vec3mer # fine, cleaned as empty list()
for (i in 1:1084){vec3mer[[i]] <- as.numeric(extrDNAkmer(as.character(mutants[[i]]), k=3, upto=T))} ; length(vec3mer); dim(vec3mer) # odel length() outputed 1084, dim() outputed NULL!

head(vec3mer) # hehe good each vector has 84 numbers, 6 made shown 
head(vec3mer, 2) # hehe good only showed 2 vectors. Each vector has 84 numbers, 
head(vec3mer, 10) # hehe good 10 vectors were shown
vec3mer # showed all
vec3mer[1] # showed [[1]]'s all, i.e., the 1st vectors of 84 numbers!


# 1.9.3 increase diversity

# 1.9.4 di-?tri-nuclueotide auto-cross covariation?xxxx


# 2.0 primary training. No need to split 6:2:2 samples first, use 10-fold cross validation in the primary training to check if selected features are good!

library(e1071) 
prime<-svm(x=vec3mer , y=df[,2] , cross=10); summary(prime) # error shit! Error in FUN(newX[, i], ...) : is.atomic(x) is not TRUE
prime<-svm(x=vec3mer[] , y=df[,2] , cross=10); summary(prime) # error shit! Error in FUN(newX[, i], ...) : is.atomic(x) is not TRUE
prime<-svm(x=vec3mer[,] , y=df[,2] , cross=10); summary(prime) # error shit! 量度数目不对
prime<-svm(x=vec3mer[[,]] , y=df[,2] , cross=10); summary(prime) # error shit! 下标数目不对
prime<-svm(x=vec3mer[[]] , y=df[,2] , cross=10); summary(prime) # error shit! 'data'的种类必需为矢量，但现在是'NULL
prime<-svm(x=vec3mer[[1:1084]] , y=df[,2] , cross=10); summary(prime) # error shit! 'Error in vec3mer[[1:1084]] : 递回索引在2层失败
prime<-svm(x=vec3mer[1:1084] , y=df[,2] , cross=10); summary(prime) # error shit! ' FUN(newX[, i], ...) : is.atomic(x) is not TRUE
prime<-svm(x=as.numeric(vec3mer) , y=df[,2] , cross=10); summary(prime) # error shit! ' 'list' object cannot be coerced to type 'double'
prime<-svm(x=as.matrix(vec3mer[[]]) , y=df[,2] , cross=10); summary(prime) # error shit! ' 'data'的种类必需为矢量，但现在是'NULL'
prime<-svm(x=as.matrix(vec3mer[]) , y=df[,2] , cross=10); summary(prime) # error shit! 'Error in FUN(newX[, i], ...) : is.atomic(x) is not TRUE

prime<-svm(x=vec3mer , y=df[,2] , cross=10); summary(prime) # error shit! 'Error in FUN(newX[, i], ...) : is.atomic(x) is not TRUE
prime<-svm(x=as.dataframe(vec3mer) , y=df[,2] , cross=10); summary(prime) # error shit! Error in as.dataframe(vec3mer) : 没有"as.dataframe"这个函
prime<-svm(x=as.data.frame(vec3mer) , y=df[,2] , cross=10); summary(prime) # error shit!  Arguments imply differing number of rows: 1084, 84
prime<-svm(x=as.vector(vec3mer) , y=df[,2] , cross=10); summary(prime) # error shit! rror in FUN(newX[, i], ...) : is.atomic(x) is not TRUE================as.vector() still useless! shit!





# fuck shit convert format first?

vecMatr<-as.matrix(vec3mer);dim(vecMatr) # odel 1084 x 1 !?
vecMatr # wrong values! shit assignment above failed!
vecMatr<-as.matrix(vec3mer[[]], nrow=1084, ncol=84);dim(vecMatr) # odel error: data'的种类必需为矢量，但现在是'NULL'
vecMatr<-as.matrix(vec3mer[], nrow=1084, ncol=84);dim(vecMatr) # odel wrong 1084 x1
vecMatr<-as.matrix(vec3mer, nrow=1084, ncol=84);dim(vecMatr) # odel wrong 1084 x1
vecMatr<-as.matrix(vec3mer[,1:84], nrow=1084, ncol=84);dim(vecMatr) # odel Error in vec3mer[, 1:84] : 量度数目不对
vecMatr<-as.matrix(vec3mer[[]][1:84], nrow=1084, ncol=84);dim(vecMatr) # odel Error in 'data'的种类必需为矢量，但现在是'NULL

vecMatr<-as.vector(vec3mer) # heyyo! no error prompted
class(vecMatr) # odel list still!
str(vecMatr) # multiple lines of num [1:84] 35282 28355 27187 35145 11539 ...
summary(vecMatr) # multiple lines of ....
length(vecMatr) # 1084
dim(vecMatr) # NULL odel!
vecMatr[1,3] # odel error !Error in vecMatr[1, 3] : 量度数目不对
vecMatr[1] # shit, seemed no diff with vec3mer list!

# online example code for converting list into vector
# copy an var first
v3<-vec3mer
v<-as.vector(unlist(v3))
class(v); str(v); # numeric type, and num ...
summary(v)
head(v) # odel only 6 were shown
length(v) # 91056!

head(v, 10) # shit, so v became an flat and continuous array/number sequences? no use for svm() either?


# below codes cannot be run! caused Rstudio into death!
prime<-svm(x=v , y=df[,2] , cross=10); summary(prime) # odel !!! out of expectation!! can run without error message shown!=========fucking stupid e1017 design! vectors "v" were not divivded into each line for each sample/instance, can it recognize each?===============bullshit 9man+ vectors, more than 30 minutes still running!===========not 30min but one-hour fuck shit non-stop!


"""heheh I was clever/smart! bullshit Rstudio running the bullshit svm()! though mouse and cursor, keys can input. But runnning failed to be terminated and Rstudio was failed to be closed!! So IT EQUALS to death of the Rstudio, failing to save everything including latest scripts and working data/variables.  My wise behavior of exporting data to xxx.RData, and script file pasted to notepad++ was very smart and correct!!!!!!  """


# reload my data saved!
load("D:/RworkDir/reconstructSVM6may2020.RData")
vec3mer # fuck ! why no data of it! not saved?
d3 # existed
seq_op # existed still
seq_ts # fuck these shits were saved, vec3mer hasn't!


# since vec3mer was not existed after restart of Rstudio, recreate it?  make the List() to be matrix() directly?

# fuck you , figured out what's wrong!  previous code's x=data, data is ok for [x,y] but now the vec3mer[x,y] will cause error! i.e., [x,y] is inaccessible for vec3mer=============yes, now convert list() type to matrix() can train normally now!



# odel libraries were missing, need to re-load "biomedR" lib for extrDNAkmer()
require(BioMedR) # reload lib
v3mat<-matrix(NA,nrow = 1084, ncol = 84)
for (i in 1:1084){v3mat[i,] <- as.numeric(extrDNAkmer(as.character(mutants[[i]]), k=3, upto=T))} ; dim(v3mat) # cool! 1084 x 84 outputed

# test
class(v3mat) # matrix da!
str(v3mat) # matrix da!
head(v3mat) # cool!
v3mat[1] # heiyo 35281
v3mat[1,1] # heiyo 35281 the same above
v3mat[1,2] # heiyo 28356
v3mat[902,2] # heiyo 28356, cool seemingly can go for svm()


# backup data/variables
save.image("D:/RworkDir/reconstructSVM6may2020.RData")

svm() # yahali, no such function, hehe recall required!

# re-call the lib
require(e1071)

prime<-svm(x=v3mat , y=df[,2] , cross=10); summary(prime) # Just used one-second ! results came out! but very bad!hehe, as shown below:
"""
Call:
svm.default(x = v3mat, y = df[, 2], cross = 10)


Parameters:
   SVM-Type:  C-classification 
 SVM-Kernel:  radial 
       cost:  1 

Number of Support Vectors:  998

 ( 499 499 )


Number of Classes:  2 

Levels: 
 benign pathogenic

10-fold cross-validation on training data:

Total Accuracy: 53.96679 
Single Accuracies:
 52.77778 48.14815 54.12844 50.92593 54.12844 62.03704 63.88889 54.12844 50.92593 48.62385 ==========bullshit accuracy!

"""


# try to shorten the features! 
prime<-svm(x=v3mat[,1:40] , y=df[,2] , cross=10); summary(prime) # odel, still very bad results! as seen blow:
"""
10-fold cross-validation on training data:

Total Accuracy: 53.96679 
Single Accuracies:
 49.07407 56.48148 55.9633 56.48148 55.9633 50 57.40741 52.29358 54.62963 51.37615 

"""




# shorten more features!
prime<-svm(x=v3mat[,1:20] , y=df[,2] , cross=10); summary(prime) #  fuck shit! no diff with above!

"""
Total Accuracy: 53.96679 
Single Accuracies:
 53.7037 53.7037 56.88073 57.40741 48.62385 50.92593 51.85185 54.12844 54.62963 57.79817 
"""
# only use first 10 features!
prime<-svm(x=v3mat[,1:10] , y=df[,2] , cross=10); summary(prime) #  fuck shit! no diff with above!

# only 5!
prime<-svm(x=v3mat[,1:5] , y=df[,2] , cross=10); summary(prime) #  fuck shit! no diff with above!

# only one-mer's 4 features!
prime<-svm(x=v3mat[,1:4] , y=df[,2] , cross=10); summary(prime) #  fuck shit! no diff with above!=============fuck masaka, primary trainin with 10-fold cross validation cannot pick the feature out?


# feature long-short switch + kernerls?
prime<-svm(x=v3mat[,1:4] , y=df[,2] , cross=10, scale=T); summary(prime) # plused scale=TRUE, shit no diff, still overall 53.96679% in average accuracy

prime<-svm(x=v3mat , y=df[,2] , cross=10, scale=T); summary(prime) # plused scale=TRUE, shit no diff, still overall 53.96679% in average accuracy



# work or argument "type="
prime<-svm(x=v3mat , y=df[,2] , cross=10, scale=T, type="nu-classification"); summary(prime) # , shit only slight diff, still overall 56.54982% in average accuracy

prime<-svm(x=v3mat , y=df[,2] , cross=10, scale=T, type="one-classification"); summary(prime) # odel worse! 50.1%!

prime<-svm(x=v3mat , y=df[,2] , cross=10, scale=T, type="eps-regression"); summary(prime) # odel errror ! said need numeric dependent variable for regression ! hehe!

# play argument "kernel="
prime<-svm(x=v3mat , y=df[,2] , cross=10, scale=T, kernel = "linear"); summary(prime) # odel 53.996679 % in average accuracy

prime<-svm(x=v3mat , y=df[,2] , cross=10, scale=T, kernel = "polynomial"); summary(prime) # odel same above 53.996679 % in average accuracy


prime<-svm(x=v3mat , y=df[,2] , cross=10, scale=T, kernel = "radial basis"); summary(prime) # odel error: wrong kernel name, hehe i typed as document indicated, your fucking stupid son!


prime<-svm(x=v3mat , y=df[,2] , cross=10, scale=T, kernel = "radial"); summary(prime) # odel still low 53.96679 %!

prime<-svm(x=v3mat[,1:4] , y=df[,2] , cross=10, scale=T, kernel = "radial"); summary(prime) # odel still low 53.96679 %!


prime<-svm(x=v3mat[,1:4] , y=df[,2] , cross=10, scale=T, kernel = "sigmoid"); summary(prime) # odel worse lower 52.8 % in average accuracy!!

tuneout<-tune(prime,train.x = v3mat , train.y = df[,2] ) # fuck shit error!
tuneout<-tune(prime,train.x = v3mat , train.y = as.character(df[,2]) ) # fuck shit still error!
tuneout<-tune(svm,train.x = v3mat , train.y = as.character(df[,2]) ) # fuck shit still error!
tuneout<-tune(svm,train.x = v3mat , train.y = df[,2] ) # odel!!!!!!!!! no errro! fuck!!!!!!!!!!!!!! the first argument tune(argu_,...) should input "svm", other shit names, variable names of primarily trained machine are not acceptable by this stupid!?
summary(tuneout) # fuck really worked as said above!!============odel but no arruracy outputed !!! BUt "Error estimation of ‘svm’ using 10-fold cross validation: 0.4602786", what your fuck is 0.46?????

# further tests!
summary(tune(svm, train.x = v3mat[,1:5], train.y = df[,2])) # still 0.4603 ! shit!
summary(tune(svm, train.x = v3mat[,1:5], train.y = df[,2]), ) # still 0.4603 ! shit! ==============if this can output 0.46??? then no need to use svm() for primary training?
summary(tune(svm, train.x = v3mat[1:500,1:5], train.y = df[1:500,2]), validation.x = v3mat[501:1000,1:5], validation.y=df[501:1000,2]) # odel error. train.x' s sample number should be greater than validation.x's?=========model is empty? fuck you? need input primarily trained svm var name?
summary(tune(prime, train.x = v3mat[1:500,1:5], train.y = df[1:500,2]), validation.x = v3mat[501:1000,1:5], validation.y=df[501:1000,2]) # odel error. 
summary(tune("prime", train.x = v3mat[1:500,1:5], train.y = df[1:500,2]), validation.x = v3mat[501:1000,1:5], validation.y=df[501:1000,2]) # odel error. "没有prime这个函数呵呵"

# directly test, but just output labels of prediction in e1071, area under curve values other performance indicators still need to call ROCR , etc 's library!

e1071::predict.svm() # error 

# caonima! tune(svm_?) tune.svm(), best.svm() differs!!?

pripred<-predict(prime, v3mat); summary(pripred) # odel feature number wrong causing error!
pripred<-predict(prime, v3mat[,1:4]); summary(pripred) # odel super bad!! yahali, 1072 benign, 12 pathogenic!hehe, sad shit!
str(pripred) # useless inf shown
pripred@ # no match 
pripred$ # no match odel then how to check the full results?

pripred # shit, error cannot show!
head(pripred) # hehe can see display!
head(pripred,1084) # max showed 1000 entries, last 84 cannot be shown!
tail(pripred,84) # cool o, i'm smart using tail(x, N-position) displyed last 84 cannot be shown!
  
10^(-6:-1) # outputed  :  1e-06 1e-05 1e-04 1e-03 1e-02 1e-01

plot(prime,data=v3mat) # error, hehe missing formula, hehe, 

tune.out<-tune(prime, 
               train.x=v3mat[,1:4], 
               train.y = df[,2], 
               validation.x = v3mat, # optional validation set?
               validation.y = df[,2], 
               ranges = list(cost = c(0.001, 0.01, 0.1, 1, 10, 100, 1000), 
               gamma = c(0.0005,0.001,0.1,0.5, 1, 2, 3, 4))); # odel  errror!  'what' must be a function or character string. Fuck trainx.=v3mat[,1:4] still error! ==========filled validation.x.y still error! shit!
summary(tune.out)



# 2.0.1 split sample 6:2:2 SET.SEED (1) is important! ========= 1 to 585 instance were benign ones; while 586 to 1084 were pathogenic! so 585 benign Versus 499 pathogenic!
# 2.0.2 firstly increase 499 pathogenic instances to 585 (i.e., 86 more are needed!), making a balanced dataset! 585*2=1170; 1170 *3/5 = 702. So in training set, 351 instances are required for each posive and negative training sets! For the rest cross validation set and testing set, 117 instances are required for each postive and negative sets!

# 不要补先？ 499patho instances分了 6:2:2先，然后只在training 过采样？防止重复的存在training 又在testing？
set.seed(1) # set.seed() has to be run EVERY TIME! NOT just once!

sample(10,20) # odel error! 不能取比总体要大的样本!

# then change splitting strategies to , 585 instances of benign/negative , splitted to 351 : 117: 117 as training:validation:testing

# postive/ patho set of 499 should be splitted to 299:100:100
# get training tags first: randomize 1 to 585 benign numbers, former 1th to 351th are as training set of benign instances; middle 352th(nd) to 468th are as cross validation set of benign instances; final 469th to 585th are as the testingset of benign instances
set.seed(2)
benign_randnumb<-sample(585,585); benign_randnumb; # cool, randomized 585 numbers! then as planned, split former, middle, and final parts!
tr_benign<-benign_randnumb[1:351]; tr_benign # as stretegized former 351 random numbers/instances are as training set of benign instances  
vali_benign<-benign_randnumb[352:468]; vali_benign # as stretegized middle 117 random numbers/instances are as cross validation set of benign instances  
test_benign<-benign_randnumb[469:585]; test_benign; length(test_benign) # as stretegized FINAL 117 random numbers/instances are as testing set of benign instances  


# similar above operation for 499 patho instances. within 1084 numbers, from 586 to 1084, 499 in total. Should be splitted to 299:100:100  
sample(10:20, 10)  # cool! worked, can specify range at the beginning, BUT! notice 10 to 20  included 11 eleven elements!!! so, size = 10, will make one missing!  

# shit here forgotten to set.seed(), though results were saved into variables! So, above set.seed(2) did not effect here!! 'coz compared sample() gave different output!
pathoRand<-sample(586:1084, 499);length(pathoRand) # cool. all elememts included! non-redundant!
# the same as benign tags' numbers' operation! 3 parts as training, validating, testing sett repectively.
trPatho<-pathoRand[1:299]; trPatho; length(trPatho) # var <- head(xx, 299 ) has the same effect?
valiPatho<-pathoRand[300:399]; valiPatho; length(valiPatho) # 100 in length, hehe, cool?
testPatho<-pathoRand[400:499]; testPatho; length(testPatho) # cool, 100.
class(pathoRand); class(test_benign) # both "integer" type var

# try to combine numbers for oversampling!
cb<-c(1,4,1,2); cb<-cbind(pathoRand[1:3],test_benign[2:4]); length(cb) # 6 in length;
cb # odel became 3 rows x 2 columns stupid! should use rbind()?
cb<-c(1,4,1,2); cb<-rbind(pathoRand[1:3],test_benign[2:4]); length(cb); cb # odel, same above, 2 dimension! shit!  6 in length;
# other simple method
cb<-c(1,35,4)
cb[3] # 4
cb[4] # NA
cb[4]<-5; cb; # cool, added 5
cb[4] # 5
cb[5]<-cb[4]; cb[5] # cool!
cb<-c(pathoRand[1:4], tr_benign[3:9]); length(cb); cb # cool!!!!!!!!! 11 in length! combined as one-dimension array!
class(cb) # "integer" cool!

""" save in the middle """
save.image("D:/RworkDir/reconstructSVM11may2020.RData")

 
# cool, just oversample by adding existed numbers . for example: tr_benign[new1_position: new2_position_number]<- tr_benign[existed positions]
# need 


# oversample patho instances, increase trainingset from 299 to 351, i.e., need to repeat 52 more, and then merge with number tags of trainingset of benign!

length(trPatho) # 299
299+52 # 351, extra 52 are required

head(trPatho,52) # fine, they are randomized numbers, add these to the end of another copy var
trpatho351<-trPatho; length(trpatho351) # 299, copied
head(trpatho351)
# addition!
trpatho351[300:351]<-head(trPatho,52); length(trpatho351); trpatho351[299:351] # odel error! shit var name wrong input,========now fixed var name fine now!
head(trPatho,52) # cool, confirmed, the same with the first 52 elements!

# integrate oversampled 351 number-tag numbers of patho instances and that of 351 benign instances
# copy the first 351 number tags of oversampled patho instances to the integrative var!
tr2x351<-trpatho351; length(tr2x351); tr2x351 # succeeded
# go on paste the next 351 instances of number-tags of benign instances' to the integrative var tr2x351
tr2x351[352:702]<-tr_benign ;  length(tr2x351) # cool, 702 in length
head(tr_benign,10); tr2x351[352:360] # cool, confirmed, the same!

tr2x351[695:702]; tail(tr_benign,10) # cool, confirmed the same!


# start primary training with ======no need? hehe tune() can directly?
# full featureset!
prm2<-svm(x=v3mat[tr2x351,] , y=df[tr2x351,2] , cross=10); summary(prm2) # shit 49.29 % in average!

prm3<-svm(x=v3mat[tr2x351,1:4] , y=df[tr2x351,2] , cross=10); summary(prm2) # shit 51.28 % in average!

# quickly insert and built an svm without cross=10, see what differs?
pr4<-svm(x=v3mat[tr2x351,], y=df[tr2x351,2]); summary(pr4) # no performance was shown
rs<-predict(pr4,v3mat[tst,]);summary(rs) # odel the same shits!!! 2 benign and 215 patho , extremely biased! failed!


# direct test with testing set, without optimization! ============= NO need to balance the positive/negative instances of testing set!? Then now they are patho:benign = 100:117, close,hehe,
# merge patho, benign two testing set's number tags first into a var
tst<-test_benign; length(tst); head(test_benign);head(tst) # no  117 element copied to the var tst. No need to worry about the order of patho and benign' coz label[position_number] will adjust order accordingly.

tst[118:217]<-testPatho; length(tst) # cool, added to 217 elements
head(tst,-1) # showed all elements except the final one whose value is 1030
tail(testPatho) # 
tst[217] # cool, 1030 matched with above!
predpr2<-predict(prm2,v3mat[tst,]); summary(predpr2) # odel! 2 benign! 215 patho! strong bias!! stupid!
predpr3<-predict(prm3,v3mat[tst,1:4]); summary(predpr3) # odel! fuck you same above horrible !! 2 benign! 215 patho! strong bias!! stupid!

""" save the data in the middle in case or loss!"""
save.image("D:/RworkDir/reconstructSVM11may2020.RData") # odel Rstudio's win path is the same with linux's / mark direction, strange!



# ==========fuck shit really tired for searching for correct usage of your fuck shit tune() functions and arguments!!=============fuck shit online example, all use the same dataset for training and cross-validation!

# integrate validation set's number-tags first also, benign : patho= 117:100, no need to balance too!?
valiset<-vali_benign;length(vali_benign);length(valiset) # cool 117 both
valiset[118:217]<-valiPatho; length(valiset) # 217, cool
tail(valiset); tail(valiPatho) # cool ,the same!



10^(-6:-1) # outputed  :  1e-06 1e-05 1e-04 1e-03 1e-02 1e-01, i.e., 0.00_xxx_1



# start tuning using full 84-length vectors' featureset
param2 <- tune(prm2, 
           train.x = v3mat[tr2x351,],
           train.y = df[tr2x351,2], 
           validation.x = v3mat[valiset,], 
           validation.y = df[valiset,2], 
           range = list(
             gamma = c(0.00001,0.00001, 0.001, 0.01,0.03,0.1,0.3,1,3,10,30),
             cost = c(0.00001,0.00001, 0.001, 0.01,0.03,0.1,0.3,1,3,10,30)), 
             tunecontrol = tune.control(sampling = "fix", fix = 1))  # fuck error! what' must be a function or character string



param2 <- tune("prm2", # changed to var name of the svm model trained, but added " " as string! 
               train.x = v3mat[tr2x351,],
               train.y = df[tr2x351,2], 
               validation.x = v3mat[valiset,], 
               validation.y = df[valiset,2], 
               range = list(
                 gamma = c(0.00001,0.00001, 0.001, 0.01,0.03,0.1,0.3,1,3,10,30),
                 cost = c(0.00001,0.00001, 0.001, 0.01,0.03,0.1,0.3,1,3,10,30)), 
               tunecontrol = tune.control(sampling = "fix", fix = 1)) #=====================fuck you! not work! no function named "prm2", hehehehe


# fuckshit, then tune(svm,___) only?
param2 <- tune(svm, # fuck really only input "svm" then no error!
               train.x = v3mat[tr2x351,],
               train.y = df[tr2x351,2], 
               validation.x = v3mat[valiset,], 
               validation.y = df[valiset,2], 
               range = list(
                 gamma = c(0.00001,0.00001, 0.001, 0.01,0.03,0.1,0.3,1,3,10,30),
                 cost = c(0.00001,0.00001, 0.001, 0.01,0.03,0.1,0.3,1,3,10,30)), 
               tunecontrol = tune.control(sampling = "fix", fix = 1)) # "tuning control = argument" is important !? Cool. finally this time without error!
summary(param2) # odel best parameters: gamma30, cost30
param2$best.parameters
param2$best.performance # 0.2396? what is this?
param2$best.model# it said this is the best model. hehe.
param2$method # "svm"
param2$nparcomb # "121"
param2$train.ind # many numbers, the randomized number tags? 
param2$sampling #  "fixed training/validation set" === what ?
param2$performances # hehe, the tuning result lists of cost and gamma combo


# fine anyway, use the best test testing set now!

pd2<-predict(param2$best.model,v3mat[tst,]); summary(pd2) ### heiyo!!!! 112 benign and 105 patho!!! suprise existed maybe!!! orginal tst set of number tags were 100 patho :117 benign!! very close o!!============further ROCR plot and computation needed!!



# start tuning, quickly test the 1:4 only 4 features!
param3 <- tune(svm, # fuck really only input "svm" then no error!
               train.x = v3mat[tr2x351,1:4],
               train.y = df[tr2x351,2], 
               validation.x = v3mat[valiset,1:4], 
               validation.y = df[valiset,2], 
               range = list(
                 gamma = c(0.00001,0.00001, 0.001, 0.01,0.03,0.1,0.3,1,3,10,30),
                 cost = c(0.00001,0.00001, 0.001, 0.01,0.03,0.1,0.3,1,3,10,30)), 
               tunecontrol = tune.control(sampling = "fix", fix = 1)) # "tuning control = argument" is important !? Cool! NO error too! 
summary(param3)  # best parameter changed differed with above o!! best gamma is 0.001,  best cost is 1e-05 ,i.e., 0.00001

# let go on predict using the testing set anyway!
pd3<- predict(param3$best.model, v3mat[tst,1:4]); summary(pd3) # odel bullshit strongly biased resuls. Still 2 benign only !! 1:4  only 4 single 1-mer feature set so shits?


# try "boot" in tune control!

param4 <- tune(svm, # fuck really only input "svm" then no error!
               train.x = v3mat[tr2x351,1:4],
               train.y = df[tr2x351,2], 
               validation.x = v3mat[valiset,1:4], 
               validation.y = df[valiset,2], 
               range = list(
                 gamma = c(0.00001,0.00001, 0.001, 0.01,0.03,0.1,0.3,1,3,10,30),
                 cost = c(0.00001,0.00001, 0.001, 0.01,0.03,0.1,0.3,1,3,10,30)), 
               tunecontrol = tune.control(sampling = "boot")) # "tuning control = argument" is important !?  NO error here, but shit "boot" type takes much longer time than "fix"! so stupid!!
summary(param4)  # best parameter changed differed with above! best gamma is 0.01 , best cost is 0.3
pd4<- predict(param4$best.model, v3mat[tst,1:4]); summary(pd4) # fuck, hehe, 2 benign, too. seriously , severely biased!===========param2 is the only hope now!




# lib ROCR compute the performance of param2
require(ROCR)
class(pd2)# factor desu,
class(df[tst,2]) # factor, too
m1<-prediction(pd2,df[tst,2]) # error in format, hehe, both formats are factor still not work, shit pack of ROCR! so stupid not smart any bit!
# still need to convert to numeric()
head(pd2)
as.numeric(pd2) # cool became 1 (benign) and 2 (patho)
head(df[tst,2]) # benign beginning
as.numeric(df[tst,2]) # cool, same above, 1 as benign

pr<-prediction(as.numeric(pd2),as.numeric(df[tst,2])) # cool. no error, to remember always need  as.numeric() hhe


pl<-performance(pr, "tpr","fpr"); plot(pl) #  sad shit!! the diagonal!!  worst performance ya!
pl@y.values # shit!!!! area under curve = 0.75!! too low ya!=======================shit, try other kernel? other classifier,e.g., randomForest? other featureset combo??===========or shuffling sets's splitting again? hehe, by luck? =========no? wrong? this is not the AUC, but stupid shit value from "tpr" and "fpr", AUC value should be "auc" in performance () 


summary(pl) # shit nothing
str(pl)


"""
1 关于svm的C以及核函数参数设置

    C一般可以选择为：10^t , t=[- 4，4]就是0.0001 到10000。选择的越大，表示对错误例惩罚程度越大，可能会导致模型过拟合================C means cost? Usually ranges from 1e0-4 to 10^4? larger value could cost overfit?? mine is 30, very large? shit?==============I have independent validation and testing set, can overcame overfitting issue? hehe.
    0）线性核函数
    （无其他参数）
    1）多项式核函数
    （重点是阶数的选择，即d，一般选择1-11：1 3 5 7 9 11，也可以选择2,4，6…）
    2）RBF核函数
    （径向基RBF内核，exp{-|xi-xj|^2/均方差}，其中均方差反映了数据波动的大小。gamma参数通常可选择下面几个数的倒数：0.1 0.2 0.4 0.6 0.8 1.6 3.2 6.4 12.8，默认的是类别数的倒数，即1/k，2分类的话就是0.5）
    3）sigmoid核函数 又叫做S形内核
    两个参数g以及r：g一般可选1 2 3 4，r选0.2 0.4 0.6 0.8 1
    4）自定义核函数



2 关于cost和gamma

SVM模型有两个非常重要的参数C与gamma。

    其中 C是惩罚系数，即对误差的宽容度。c越高，说明越不能容忍出现误差,容易过拟合。C越小，容易欠拟合。C过大或过小，泛化能力变差

    gamma是选择RBF函数作为kernel后，该函数自带的一个参数。隐含地决定了数据映射到新的特征空间后的分布，gamma越大，支持向量越少，gamma值越小，支持向量越多。支持向量的个数影响训练与预测的速度。
    Grid Search
    使用grid Search虽然比较简单，而且看起来很naïve。但是他确实有两个优点：
    可以得到全局最优
    (C,gamma)相互独立，便于并行化进行

"""

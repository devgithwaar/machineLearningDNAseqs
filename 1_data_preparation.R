
# R codes

# Author: Hui-Heng Lin, PhD. Copyright reserved. 12th feb. 2021

# Purpose: DNA sequence Data preparation---- to prepare /generate BRCA1 gene variants/mutants for machine learning analysis



# Notice that because gene mutation/variant databases like ensembl and NCBI-clinVar DO NOT PROVIDE THE FULL-LENGTH mutated gene DNA sequences, instead, they only display and provide mutation/variation data such as the mutation position and seqeuence changes --- i.e., we have to download the reference/standard full-length gene sequence and then manually create the mutated full-lengt DNA sequences via coding.


""" 

BRCA1 gene mutants' preparation

1. Download raw/standard BRCA1 full length gene(DNA sequences) from genome browser (Ensembl, NCBI, etc) displaying the human reference genome GCh38.


2. Remember and save the start coordinate of BRCA1 gene in Gch38 version of human genome, which is 43044285


3. Save the BRCA1 full length gene DNA sequences to local computer .txt file (or alternatively certain genome browser may provide users with .fasta format file) and loaded it into R computing environment.

4. Query ensembl or NCBI-clinVar database for BRCA1 gene's mutation information. --- in this work, we only handle cases of single point/base mutation of BRCA1 gene, other types of more complexed mutation forms are not considered here.

5. Collect the mutation information including the mutation positions, pathogenic risks(benign or pathogenic) and the base substitutions, and then store them in the file. (may be in .csv format)

"""
 
s1<-43044285; print(s1); print(s1-3) # The start coordinate of BRCA1 gene in Gch38 version of human genome, which is 43044285. Assign to variable for later mutation operations


df<-read.csv("I:/your_path_way/mutation_info.csv", header=T)  # from  disk file load the numbers, BRCA1 mutations/SNPs location/genome co-ordinates, mutatated bases, and label(benign or patho) into R environment

df;  # Check the data. loaded with number, mutation location/co-ordinates, mutatated bases, label:benign/patho

str(df) # the data_frame type. Mutated bases should be df[,4]; co-ordinates should be df[,3], label of pathogenic/benign are df[,2] 



#  test access of data via [x,y]
df[1,3] # the co-ordinate was , 43094550

df[2,4] # the nucleotide base was the "C"

df[7,3] # 43071102，the co-ordinate




library(Biostrings) # need to use the DNAString() for format conversion

dna<-scan( file="I:/your_path_/fulllengthBRCA1_DNA_sequence.txt", what = character(0) ) # 
# load the full sequence into R environment via basic function scan( ) . Note that the raw format of DNA sequence might contain  multiple "\n" new line marks. It is necessary to remove them and make the full length sequence into one-line format.

head(dna) ; class(dna) # checking the DNA sequence imported. The variable was in character data type

dna[1] # access the header annoation of the sequence:">GCh38 human reference genome..."

dna[2] # acces the full length BRCA1 gene sequence. Alternatively, dna[2][1] can access the same content as well.

dna[3] # NA


# Convert the sequence to "biostrings" data type
DNAString(dna[2]) # check

d3<-DNAString(dna[2]) # assign the converted format to variable

class(d3) # data type checking. Has become the "biostrings" type. Conversion succeeded.

d3[2]; d3[1];d3[2]; d3[4]  # through indication of seqeuence position number, we could access any single base of the full length sequence. Prepare for the mutation operation (base substitution). 

substr(d3,1,2) #  way to access multiple neighbor DNA bases


#  copy another variable for trials of DNA single base mutation operation
seq_ts<-d3; 

# check the data 
class(seq_ts) # data type: "biostrings"


# confirm different ways for  accessing DNA single base 
seq_ts[5] 

seq_ts[43071102-s1] # s1 is aforementioned start co-ordinate of BRCA1 gene in the genome

seq_ts[df[7,3]-s1] 


# mutation operation trial
seq_ts[5]<-"K"; seq_ts[5] # base was substituted to K 

seq_ts[5]<-"T"; seq_ts[5] # reverse the base to the original one "T" 



"""
For creating mutants and the final numeric vectors/representations. Two options are available.

First one is to generate all the mutant seqs first, and then convert them into vectors afterwards; 

The second one is using a loop that creating one mutant sequence, and then sebsequently generate its numeric representation right away and store into a variable. 

Option 1 is the less complex way to go.

"""




# generate the same type: "biostrings" for mutation operations 
seq_op<-d3; 







mutants<-list() # firstly, create an emtpy list var is required to store each mutant seq


# the loop for creating mutants

for (i in 1:y) { # y is your number of mutants
  
    seq_op[df[i,3]-s1+1]<-as.character(df[i,4]); # replace one base each time, Note that "as.character( )" is required here. Otherwise error.  s1 is aforementioned start co-ordinate of BRCA1 gene in the genome, hence the number of base mutation position -s1 is the number of actual mutation position of BRCA1 gene (rather than the full genome's )   
    
    mutants[i] <- seq_op; # save the mutant into the list varible

    seq_op <- d3 # every time 's ending, recover it to wildtype/original standard BRCA1 seq is required!
};  






# check the output result of above loop, so as to confirm if mutant seqs were generated successfully

dim(mutants) #  NULL because it is a list tyep

length(mutants) # output the length of the list variable

mutants # list output

mutants[1] #  the 1st mutated seq was shown

mutants[[1]] # the same above




# check if mutations were successful

df[7,4] # check the original mutated base from database

substr(mutants[[7]],43071102-s1,43071102-s1) # locate the mutation spot and mutated base, and compare with above to see if they are consistent.




#  save/backup in the middle: all var into .Rdata. Backup is vital!  can either do it via clicking graphical interface on Rstudio or use command:

save.image(file = "D:/your_path_folder/mutant_seq_generation.RData", version = 2.0 , safe = TRUE) 


import os
from Bio import SeqIO, Entrez

#os.mkdir("MiniProject_Delaney_Sauer")
os.chdir("MiniProject_Delaney_Sauer")
outfile = open("miniProjectLog.log",'w')

##1. First, retrieve the following transcriptomes from two patient donors
##from SRA and convert to paired-end fastq files. You can use wget
##(by constructing the path based on the SRR numbers for each of these samples).

#PlAN: use os.system to call unix shell, have list of websites to use wget on
#d12="https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/SRR5660030/SRR5660030.1"
#d16="https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/SRR5660033/SRR5660033.1"
#d22="https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/SRR5660044/SRR5660044.1"
#d26="https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/SRR5660045/SRR5660045.1"
    
#os.system("wget " + d12)
#os.system("wget " + d16)
#os.system("wget " + d22)
#os.system("wget " + d26)
#files are under the last SRX number in the url 
#need to create fastq files using fastq dump
#os.system("fastq-dump -I --split-files SRR5660030.1")
#os.system("fastq-dump -I --split-files SRR5660033.1")
#os.system("fastq-dump -I --split-files SRR5660044.1")
#os.system("fastq-dump -I --split-files SRR5660045.1")
#both the wget and dump could have been done in a loop, I wanted to 
#code each line to make debugging easier
    
##2. Build a transcriptome index for HCMV (NCBI accession EF999921).
##Use Biopython to retrieve input and then build the index with kallisto.

Entrez.email = "delaneyjosauer@gmail.com.com"
term = 'EF999921'
handle = Entrez.efetch(db="nucleotide", id=[term], rettype="gb") #searching genbak via Entrez
record = SeqIO.read(handle, "genbank")#read the genbank file

count = 0
index = open("kalIndex.fasta","w")
for feature in record.features:
    if feature.type == 'CDS':
        CDS = str(feature.qualifiers['translation']).replace("'","")
        CDS = CDS.replace("[","")
        CDS = CDS.replace("]","")
        #the .replace methods are nessecary to format the CDS to a fasta file
        name = str(feature.qualifiers['protein_id']).replace("'","")
        name = name.replace("[","")
        name = name.replace("]","")
        CDS = ">" + name + "\n" + CDS #makes the seqeunces into fasta format
        index.write(CDS)
        count += 1
outfile.write("The HCMV genome (EF999921) has " + str(count) + " CDS.")

##3. Quantify the TPM of each CDS in each transcriptome using kallisto and use these results as input to find differentially expressed
##genes between the two timepoints (2pi and 6dpi) using the R package sleuth.Write the following details for each significant
##transcript (FDR < 0.05) to your log file, include a header row, and tab-delimit each item: 

#okay let's build the index
os.system("kallisto index -i index.idx kalIndex.fasta")
#does the index need to be nucleotide? its aa from genbank
os.system("kallisto quant -i index.idx -o SRR5660030 -b 30 -t 2 SRR5660030.1_1.fastq  SRR5660030.1_2.fastq") 
os.system("kallisto quant -i index.idx -o SRR5660033 -b 30 -t 2 SRR5660033.1_1.fastq  SRR5660033.1_2.fastq")
os.system("kallisto quant -i index.idx -o SRR5660044 -b 30 -t 2 SRR5660044.1_1.fastq  SRR5660044.1_2.fastq")
os.system("kallisto quant -i index.idx -o SRR5660045 -b 30 -t 2 SRR5660045.1_1.fastq  SRR5660045.1_2.fastq")

#to run Sleuth, need an prewritten r script that we're going to access from my github
#insert wget method here
#os.system("Rscript MiniProjectCode.R")
#write to log the return

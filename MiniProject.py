import os
from Bio import SeqIO, Entrez

os.mkdir("MiniProject_Delaney_Sauer")
#os.chdir("MiniProject_Delaney_Sauer")
outfile = open("miniProjectLog.log",'w')

##1. First, retrieve the following transcriptomes from two patient donors
##from SRA and convert to paired-end fastq files. You can use wget
##(by constructing the path based on the SRR numbers for each of these samples).

#PlAN: use os.system to call unix shell, have list of websites to use wget on
d12="https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/SRR5660030/SRR5660030.1"
d16="https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/SRR5660033/SRR5660033.1"
d22="https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/SRR5660044/SRR5660044.1"
d26="https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/SRR5660045/SRR5660045.1"
    
os.system("wget " + d12)
os.system("wget " + d16)
os.system("wget " + d22)
os.system("wget " + d26)
 
#need to create fastq files using fastq dump
os.system("fastq-dump -I --split-files SRR5660030.1")
os.system("fastq-dump -I --split-files SRR5660033.1")
os.system("fastq-dump -I --split-files SRR5660044.1")
os.system("fastq-dump -I --split-files SRR5660045.1")
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
		CDS = str(feature.location.extract(record).seq) #get sequence of CDS
		name = str(feature.qualifiers['protein_id']).replace("'","") #use string methods to format 
		name = name.replace("[","")
		name = name.replace("]","")
		CDS = ">" + name + "\n" + CDS #makes the seqeunces into fasta format
		index.write(CDS + "\n") #makes kalIndex fasta
		count += 1
outfile.write("The HCMV genome (EF999921) has " + str(count) + " CDS." + "\n")

##3. Quantify the TPM of each CDS in each transcriptome using kallisto and use these results as input to find differentially expressed
##genes between the two timepoints (2pi and 6dpi) using the R package sleuth.Write the following details for each significant
##transcript (FDR < 0.05) to your log file, include a header row, and tab-delimit each item: 

#okay let's build the index
os.system("kallisto index -i index.idx kalIndex.fasta")
os.system("kallisto quant -i index.idx -o SRR5660030 -b 30 -t 2 SRR5660030.1_1.fastq  SRR5660030.1_2.fastq") 
os.system("kallisto quant -i index.idx -o SRR5660033 -b 30 -t 2 SRR5660033.1_1.fastq  SRR5660033.1_2.fastq")
os.system("kallisto quant -i index.idx -o SRR5660044 -b 30 -t 2 SRR5660044.1_1.fastq  SRR5660044.1_2.fastq")
os.system("kallisto quant -i index.idx -o SRR5660045 -b 30 -t 2 SRR5660045.1_1.fastq  SRR5660045.1_2.fastq")

#need to create text file with paths for r script
kalTable = open("kalTable.txt",'w')
kalTable.write("sample    condition    path" + "\n")
kalTable.write("SRR5660030 D12 SRR5660030" + "\n")
kalTable.write("SRR5660033 D16 SRR5660033" + "\n")
kalTable.write("SRR5660044 D22 SRR5660044" + "\n")
kalTable.write("SRR5660045 D26 SRR5660045" + "\n")

#3/2 this does not work
#to run Sleuth, need an prewritten r script that we're going to access from my github
#os.system("wget https://github.com/d-sauer-9/MiniProject/blob/main/SleuthScript.R")
#os.system("Rscript SleuthScript.R")
#write to log the return
outfile.write("As of 3/4, Sleuth component does not work")

#running Bowtie 2
os.system("wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/845/245/GCF_000845245.1_ViralProj14559/GCF_000845245.1_ViralProj14559_genomic.fna.gz")
#GCF_000845245.1_ViralProj14559_genomic.fna.gz to create index
os.system("bowtie2-build GCF_000845245.1_ViralProj14559_genomic.fna.gz HCMVRef")
#do the mapping
os.system("bowtie2 --quiet -x HCMVRef -1 SRR5660030.1_1.fastq -2 SRR5660030.1_2.fastq -S Donor1Day2.sam -S --al-conc-gz SRR5660030.1_mapped%.fq.gz")
os.system("bowtie2 --quiet -x HCMVRef -1 SRR5660033.1_1.fastq -2 SRR5660033.1_2.fastq -S Donor1Day6.sam -S --al-conc-gz SRR5660033.1_mapped%.fq.gz")
os.system("bowtie2 --quiet -x HCMVRef -1 SRR5660044.1_1.fastq -2 SRR5660044.1_2.fastq -S Donor2Day2.sam -S --al-conc-gz SRR5660044.1_mapped%.fq.gz")
os.system("bowtie2 --quiet -x HCMVRef -1 SRR5660045.1_1.fastq -2 SRR5660045.1_2.fastq -S Donor2Day6.sam -S --al-conc-gz SRR5660045.1_mapped%.fq.gz")
#this takes a long amount of time

#find out to get num of reads before and after and write to log file
#use loop to find the num of reads before and after
countd12  = 0
countd122 = 0
for record in SeqIO.parse("SRR5660030.1_1.fastq", "fastq"):
	countd12+=1
outfile.write("Donor 1 (2dpi) had " + str(countd12) + " reads before Bowtie2 and " + str(countd122) +  " reads after" +"\n")

countd16  = 0
count162 = 0
for record in SeqIO.parse("SRR5660033.1_1.fastq", "fastq"):
	countd16+=1
outfile.write("Donor 1 (6dpi) had " + str(countd16) + " reads before Bowtie2 and " + str(count162) + " reads after" + "\n")

countd22  = 0
countd222 = 0
for record in SeqIO.parse("SRR5660044.1_1.fastq", "fastq"):
	countd22 += 1
outfile.write("Donor 2 (2dpi) had " + str(countd22) + " reads before Bowtie2 and " + str(countd222) + " after"  +"\n")

countd26  = 0
countd262 = 0
for record in SeqIO.parse("SRR5660045.1_1.fastq", "fastq"):
	countd26+=1
outfile.write("Donor 2 (6dpi) had " + str(countd26) + " reads before Bowtie2 and " + str(countd262) + " reads after"  +"\n")

#using Bowtie2 reads in SPAdes
os.system("spades -k 55,77,99,127 -t 2 --only-assembler -s SRR5660030.1_mapped%.fq.gz -s SRR5660033.1_mapped%.fq.gz -s SRR5660044.1_mapped%.fq.gz -s SRR5660045.1_mapped%.fq.gz -o assembly/")
outfile.write("SPAdes Command used: spades -k 55,77,99,127 -t 2 --only-assembler -s SRR5660030.1_mapped%.fq.gz -s SRR5660033.1_mapped%.fq.gz -s SRR5660044.1_mapped%.fq.gz -s SRR5660045.1_mapped%.fq.gz -o assembly/" + "\n")

#finding number of contigs > than 1000
file = "assembly/contigs.fasta"
fileread = SeqIO.parse(file,'fasta')
count = 0
total = 0
longestContig = ""
longestContigLen = 1
for entry in fileread:
	temp = len(entry.seq)
	if temp>= 1000:
		count += 1
		total = total + temp
	if temp > longestContigLen:
		longestContig = str(entry.seq)
outfile.write("There are " + str(count) + " contigs > 1000 bp in the assembly" + "\n")
outfile.write("There are " + str(total) + " bp in the assembly" +"\n")
#make longest contig a file to blast in fasta format
fastaMaker = open("longestCongtig.fasta","w")
fastaMaker.write("> Longest Contig Read" + "\n" + str(longestContig))

#using BLAST
#for some reason, makeblastdb does not like the file in my repo, so I'm going to copy it into another file
database = open("database.fasta","w")
reader = SeqIO.parse("databaseFasta.fasta",'fasta')
for entry in reader:
	database.write(">" + str(entry.id) + "\n" + str(entry.seq) + "\n")

#make local database
#as of 3/4 not working; says database.fasta is not a fasta file
#os.system("makeblastdb -in database.fasta -out database -title database -dbtype nucl")
#os.system("blastn -query longestContig.fasta -db database out myResults.csv -outfmt 7")
outfile.write("As of 3/4, BLAST component does not work")

outfile.close()
reader.close()
database.close()
fastaMaker.close()
fileread.close()
kalTable.close()

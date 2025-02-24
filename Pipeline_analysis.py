#argparse code
import argparse
import sys
import os
import subprocess
import numpy as np
import pandas as pd
import zipfile
import glob
import shutil
from Bio import SeqIO
from Bio import Entrez


#FIRST setting up the directory that all the output files should go to
output_dir = "/home/2025/aavalos4/PythonPipeline_Lexi_Avalos" #change directory here

#Creating the directory if it isn't there already
os.system(f"mkdir -p {output_dir}")

#move into the directory from now on
os.chdir(output_dir)



#Step 2: Quantify TPM in each sample using Kallisto
# a) set up and get info from NCBI
Entrez.email = "aavalos4@luc.edu" #NCBI needs email

#defining variables to store files and accession number
accession = "NC_006273.2" #the accession number for the HCMV genome
HCMV_fasta_file = os.path.join(output_dir, "HCMV_reference.fasta") #the output place where the transcriptome fasta will be put - from (2-b)
log_file = os.path.join(output_dir, "PipelineProject.log") #the log file to hold all the details throughout the project
kallisto_index = os.path.join(output_dir, "HCMV_index_kallisto") #output for the kallisto index from (2-c)

#Get the HCMV genome from NCBI
#database = nucleotide, id = NC_006273.2, rettype = Genbank format, format/retmode = text
handle = Entrez.efetch(db='nucleotide', id = accession, rettype = "gb", retmode= "text")
record = SeqIO.read(handle, "genbank") #making the downloaded Genbank record into an object for easy use
#print(record) #checking how it looks
handle.close() #closing connection

# b) Extracting the coding sequence (CDS) features from the HCMV genome - to use as input for kallisto index command
cds_count = 0 # to count how many cds sequences are taken out
with open(HCMV_fasta_file, "w") as fasta_output: #opening the HCMV_fasta_file, and in write mode so it writes to this file
    for variable in record.features: #going through all the features in the record
        if variable.type == "CDS" and "protein_id" in variable.qualifiers: #if the the feature is cds, make sure it has a protein id 
            protein_id = variable.qualifiers["protein_id"][0] #getting the protein id to put as the header CDS identifier
            cds_sequence = variable.extract(record.seq) #getting CDS sequence from record object
            fasta_output.write(f">{protein_id}\n{cds_sequence}\n") #writes the cds into the fasta file with the correct format
            cds_count += 1 #will add one to the count every time it extracts a cds sequence

#print(cds_count) #checking

#putting the cds info into the log
with open(log_file, "w") as log_output:
    log_output.write(f"The HCMV genome ({accession}) has {cds_count} CDS.\n") #making to put a newline so the next thing isn't conjoined
    log_output.write(f"\n") #putting a newline for more space


# c) Build a transcriptome index for HCMV (NCBI accession NC_006273.2)
#Making the index for HCMV
making_kallisto_index = f'kallisto index -i {kallisto_index} {HCMV_fasta_file}' #command for it to run
#does kallisto index command - shell = True so it runs good bc i wrote command as a string
subprocess.run(making_kallisto_index, shell = True)
#print(making_kallisto_index) #checking to see what it prints






#Step 3: Quantify the TPM of each CDS in each transcriptome using kallisto 
# a) get SRR files 
srr_pathway = os.path.join(output_dir, "SRR_file") #change for new Srr
os.makedirs(srr_pathway, exist_ok=True) #make sure it makes it

#moving all SRR files to the SRR_file
for file in os.listdir(output_dir): #making sure of all the files in the output directory
    source = os.path.join(output_dir, file) #starting place
    place = os.path.join(srr_pathway, file) #place we want to move it to
    if file.startswith("SRR") and "fastq" not in file: #making sure it is the right one
        if os.path.isfile(source): #making sure it is there; only moving files not directories- debug
            shutil.move(source, place) #moving it

#removes whitespace, gets a list of all the files in directory, and only chooses the files with the name that starts with SRR
srr_names = [filename.strip() for filename in os.listdir(srr_pathway) if filename.startswith("SRR")]

#print(srr_names) #checking

# b) get fastq files
fastq_pathway = os.path.join(output_dir, "fastq") #change for new fastqs
os.makedirs(fastq_pathway, exist_ok=True) #makes sure it makes it

#doing the same thing to move fastq files to their directory
for file in os.listdir(output_dir): #checking for files
    if file.endswith(".fastq"): #moving only fastq files
        source = os.path.join(output_dir, file) #starting point
        place = os.path.join(fastq_pathway, file) #ending point
        if os.path.exists(source): #mkaing sure it is there
            shutil.move(source, place) #moving it 

#making a dictionary called fastq_files where the key is srr and value is a tuple that has paths for fastq files - associating fastq with their SRR
fastq_files = {srr: (os.path.join(fastq_pathway, f"{srr}_1.fastq"), #os.path.join = makes path for first fastq file
					os.path.join(fastq_pathway, f"{srr}_2.fastq")) #same but for the second fastq file
				for srr in srr_names} #looping through all SRR numbers in SRR file

print(fastq_files) #checking

# c) run kallisto
kallisto_output = os.path.join(output_dir, "kallisto_results") #place to put kallisto output

for srr, (fastq1, fastq2) in fastq_files.items(): #going through fastq files dictionary 
	kal_output = os.path.join(kallisto_output, srr) #mkaes a separate directory for each sample in the output folder
	os.makedirs(kal_output, exist_ok= True) #makes directories and debug (exist_ok = True) will work if directory already is there

	#kallisto command - using the previous made index as input, outputting the new output, bootstraps = 5 for faster use, using both fastqs and a thread speed of 4
	command_kallisto = f"time kallisto quant -i {kallisto_index} -o {kal_output} -b 30 --threads 4 {fastq1} {fastq2}"
	subprocess.run(command_kallisto, shell = True) #running the kallisto command

# d) get the conditions from homemade txt file
conditions = os.path.join(output_dir, "conditions.txt")

srr_conditions = {} #empty dictionary to hold conditions and SRR numbers

with open(conditions, "r") as file: #open the conditions txt file
    for line in file: #for each line in the file 
        sections = line.strip().split(": ") #splitting at the colon
        if len(sections) == 2: #which it should because there are only 2 sides of the colon
            condition, srr = sections #getting condition and srr numbers
            srr_conditions[srr] = condition #making the srr number as the key and the condition as the value
            
#print(srr_conditions) #checking
            

# e) Get minimum, median, mean, and maximum TPM from the abundance.tsv kallisto output file
abundance_tsv = os.path.join(output_dir, "kallisto_results/{srr}/abundance.tsv") #pathway to get abundance_tsv file

tpm_data = [] #empty list to put all tpm data

for srr, condition in srr_conditions.items(): # for each srr number and condition in the dictionary
    abundance_tsv = abundance_tsv.format(srr = srr) #changing srr with actual srr name
    if os.path.exists(abundance_tsv): #if the path for abundance tsv exists
        read_tsv = pd.read_csv(abundance_tsv, sep="\t") #read the tsv file 
        if "tpm" in read_tsv.columns: #making sure tpm column is there
            tpm_numbers = read_tsv["tpm"].values #get tpm values
            
			#getting all the calculations from tpm_numbers
            min_tpm = np.min(tpm_numbers) #the minimum
            med_tpm = np.median(tpm_numbers) #the median
            mean_tpm = np.mean(tpm_numbers) #the mean
            max_tpm = np.max(tpm_numbers) #the max
            
			#store into tpm list
            tpm_data.append([srr, condition, min_tpm, med_tpm, mean_tpm, max_tpm])
            #print([srr], condition, min_tpm, med_tpm, mean_tpm, max_tpm)

#log_file = "/home/2025/aavalos4/PythonPipeline_Lexi_Avalos/PipelineProject.log"

#putting stats of tpm to log file
with open(log_file, "a") as output: #want to append instead of writing over it 
    output.write("sample\tcondition\tmin_tpm\tmed_tpm\tmean_tpm\tmax_tpm\n") #making the header row
    for i in tpm_data: #for each variable in tpm_data
        #changing list into strings, joins all of them by a tab and then makes a newline
        output.write("\t".join(map(str, i))+ "\n") #writing it to the file
    output.write(f"\n") #extra space for the next thing

    


#Step 4: Use output from kallisto as input for Sleuth

# a) make sleuth table that it will accept
#need a table of samples, conditions and paths to kallisto results
samples = srr_conditions #samples and their conditions were already in a dictionary from previous step
#kallisto_output already has the kallisto results
sleuth_input_file = os.path.join(output_dir, "sleuth_input.txt") #output file pathway

with open(sleuth_input_file, "w") as input:
    input.write("sample\tcondition\tpath\n") #writing the header
    for srr, condition in samples.items(): # for each srr number and condition in the samples dictionary
        path_sample = os.path.join(kallisto_output, srr) #the path to the kallisto results
        input.write(f"{srr}\t{condition}\t{path_sample}\n") #write the SRR number, condition, and sample path to the kallisto results

# b) run the sleuth analysis using the sleuth.R wrapper
rscript_sleuth = os.path.join(output_dir, "sleuth.R") #pathway to sleuth.R wrapper

sleuth_output = os.path.join(output_dir, "sleuth_output.tsv") #output pathway

#running the sleuth rscript
subprocess.run(["Rscript", rscript_sleuth], check = True) #running the rscript

#add the results to the log file
if os.path.exists(sleuth_output): #if there is a path for sleuth results
    with open(sleuth_output, "r") as result, open(log_file, "a") as output: #read sleuth output, append into the log file
        for line in result: #for each line in the sleuth output
            output.write(line) #append each of the sleuth results to the log
        output.write(f"\n") #creating another line for more space





#Step 5: Let's find the strains that are most similar to the patient samples - part 1 (assembly with Bowtie2)

# a) Function to saving the amount of reads before and after Bowtie2 filtering - easier for faster speed
def counting_reads(fastq_file): #take in fastq path to run in this function
    if not os.path.exists(fastq_file): #trying to debug - if the path doesn't exist
        print(f"fastq file {fastq_file} not found") #put fastq file as not found
        return 0 #giving back 0
    
    #wc l is counting the number of lines in the fastq; capture_output is making sure to store output; and text = output is string instead of bytes
    reads_num = subprocess.run(["wc", "-l", fastq_file], capture_output=True, text = True)
    return int(reads_num.stdout.split()[0]) // 4 #taking the first value (the num of reads), int changes into integer, and /4 bc fastq has 1 read per every 4 lines
    
# b) Using Bowtie2, creating a genome index for HCMV (NC_006273.2)
#the genome from NCBI is under the variable (record) from before
bowtie2_record = record #the hcmv reference fasta
#print(bowtie2_record) #checking

#defining new output fasta file
bowtie2_hcmv_fasta = os.path.join(output_dir, "bowtie2_hcmv_fasta")

#taking out the sequence from record to the bowtie2 hcmv fasta file
with open(bowtie2_hcmv_fasta, "w") as output:
    output.write(f">{record.id} {record.description}\n") #writing the fasta header with the id and description
    output.write(str(record.seq) + "\n") #writing the sequence

#pathway to put the bowtie2 index with index prefix 
bowtie2_index = os.path.join(output_dir, "bowtie2_index")

#run bowtie2 index
subprocess.run(["bowtie2-build", bowtie2_hcmv_fasta, bowtie2_index], check = True)


# c) Running bowtie2
#have fastq files and srr in this one 
#print(fastq_files) 

#output pathway for bowtie2 results
bowtie2_output = os.path.join(output_dir, "bowtie2_output")
os.makedirs(bowtie2_output, exist_ok=True) #makes a directory bc will fail if not already there

#for correct format in log file, need to add in the full conditions
full_conditions = os.path.join(output_dir, "full_conditions.txt") #path to txt file with full conditions

#dictionary with full conditions
dict_full_conditions = {} #empty for it to have it hold info

#reading the file and storing it
with open(full_conditions, "r") as file: #opening and reading the txt file
    for line in file: #for each line in the file
        sections = line.strip().split(": ") #splitting at the colon
        if len(sections) == 2: #if more than 2 sections, wont work
            condition, srr = sections #these two categories make up the sections
            dict_full_conditions[srr] = condition #srr as the key and conditions as the value


#adding all of this to my log file
with open(log_file, "a") as output:
    #going through fastq_files dictionary and running Bowtie2
    for srr, (fastq1, fastq2) in fastq_files.items(): #for each srr number and fastq1 and 2 in fastq_files dictionary
        sam_output = os.path.join(bowtie2_output, f"{srr}.sam") #add the sam output with this name to the bowtie2_output bc gives out a .sam file
        mapped_fastq_prefix = os.path.join(bowtie2_output, f"{srr}_mapped") #making the prefix for the fastq files that come out as bowtie2 output
        map_fastq1 = f"{mapped_fastq_prefix}_1.fq.gz" #mapped 1st read -- using this specificically for counting reads after bowtie2 run
        map_fastq2 = f"{mapped_fastq_prefix}_2.fq.gz" #mapped 2nd read
        mapped_fastq = f"{mapped_fastq_prefix}_%.fq.gz" #bowtie2 format when outputted

        #counting the reads before running the bowtie2 command

        reads_before = counting_reads(fastq1) #calling the function a lil further up north

        #print(reads_before)
        #running the bowtie2 command
        subprocess.run(["bowtie2", "--quiet", "-x", bowtie2_index, "-1", fastq1, "-2", fastq2, "-S", sam_output, "--al-conc-gz", mapped_fastq], check = True)

        #counting the reads after running the bowtie2 command
        reads_after = counting_reads(map_fastq1) #getting mapped read count 

        #getting the condition to use when doing the log
        condition = dict_full_conditions.get(srr, srr) #going to just put srr if not found

        #getting the results into the log
        output.write(f"{condition} had {reads_before} read pairs before Bowtie2 filtering and {reads_after} read pairs after. \n")
    output.write(f"\n") #creating another line for extra space


'''   
#ignore this section -- this was so i wouldnt have to keep redoing the bowtie part bc it took so long
log_file = "/home/2025/aavalos4/PythonPipeline_Lexi_Avalos/PipelineProject.log"
# a) get SRR files 
srr_pathway = "/home/2025/aavalos4/PythonPipeline_Lexi_Avalos/SRR_file" #change for new Srr
#removes whitespace, gets a list of all the files in directory, and only chooses the files with the name that starts with SRR
srr_names = [filename.strip() for filename in os.listdir(srr_pathway) if filename.startswith("SRR")]

#print(srr_names) #checking

# b) get fastq files
fastq_pathway = "/home/2025/aavalos4/PythonPipeline_Lexi_Avalos/fastq" #change for new fastqs
#making a dictionary called fastq_files where the key is srr and value is a tuple that has paths for fastq files - associating fastq with their SRR
fastq_files = {srr: (os.path.join(fastq_pathway, f"{srr}_1.fastq"), #os.path.join = makes path for first fastq file
					os.path.join(fastq_pathway, f"{srr}_2.fastq")) #same but for the second fastq file
				for srr in srr_names} #looping through all SRR numbers in SRR file

#bowtie2 output
bowtie2_output = "/home/2025/aavalos4/PythonPipeline_Lexi_Avalos/bowtie2_output"
'''




#Step 6: Using bowtie2 output reads fastqs in SPAdes to generate 2 assemblies - one for each patient/donor  - using kmer size 77
#making a pathway and place to put spades output
spades_output = os.path.join(output_dir, "spades_output")
os.makedirs(spades_output, exist_ok = True) #making this just in case it gives error

#initialize the kmer size
kmer_size = "77" #string so i dont get error

#add to log file while running
with open(log_file, "a") as output:
    #we want to go through all of the samples - so make a loop
    for srr in fastq_files.keys(): #for each srr in the dictionary
        spades_name = os.path.join(spades_output, f"{srr}_assembly") #putting the name for the output of spades

        #get mapped fastqs when not trying to run bowtie2 again
        map_fastq1 = os.path.join(bowtie2_output, f"{srr}_mapped_1.fq.gz")
        map_fastq2 = os.path.join(bowtie2_output, f"{srr}_mapped_2.fq.gz")

        #making a variable for spades command to print in log file -- this is the command in a variable
        spades_run = ["spades.py", "-k", kmer_size, "-t", "2", "--only-assembler", "-1", map_fastq1, "-2", map_fastq2, "-o", spades_name]

        #running the spades command
        subprocess.run(spades_run, check = True)

        #putting spades command into log
        output.write(f"Spades for {srr}: {' '.join(spades_run)}\n") #going to show the command for each SRR number
    output.write(f"\n") #extra line for extra space

'''
#ignore this section - spades took a while
# a) get SRR files 
srr_pathway = "/home/2025/aavalos4/PythonPipeline_Lexi_Avalos/SRR_file" #change for new Srr
#removes whitespace, gets a list of all the files in directory, and only chooses the files with the name that starts with SRR
srr_names = [filename.strip() for filename in os.listdir(srr_pathway) if filename.startswith("SRR")]

#print(srr_names) #checking

# b) get fastq files
fastq_pathway = "/home/2025/aavalos4/PythonPipeline_Lexi_Avalos/fastq" #change for new fastqs
#making a dictionary called fastq_files where the key is srr and value is a tuple that has paths for fastq files - associating fastq with their SRR
fastq_files = {srr: (os.path.join(fastq_pathway, f"{srr}_1.fastq"), #os.path.join = makes path for first fastq file
					os.path.join(fastq_pathway, f"{srr}_2.fastq")) #same but for the second fastq file
				for srr in srr_names} #looping through all SRR numbers in SRR file

#spades output
spades_output = "/home/2025/aavalos4/PythonPipeline_Lexi_Avalos/spades_output"

'''


#Step 7: Finding the strains that each assembly aligns to 

# a) Get the longest contig from each SPAdes assembly
long_contig = {} #making an empty dictionary to hold the longest contig from each spades assembly

for srr in fastq_files.keys(): #going through each SRR number in spades output
    contigs_info = os.path.join(spades_output, f"{srr}_assembly", "contigs.fasta") #the path to get to the spades output with the contigs
    
    if os.path.exists(contigs_info): #making sure contigs.fasta is there
        longest_contig = None #starting at none and it will be replaces to hold longest contig
        longest_contig_num = 0 #place so it can count to see if the contig is longer or not

    #get the contig from the file
    for contig in SeqIO.parse(contigs_info, "fasta"): #for each contig in the fasta file
        contig_length = len(contig.seq) #getting the length of the contig in use
        if contig_length > longest_contig_num: #if the contig length is longer than the other contigs
            longest_contig = contig #makign the whole object as the full seq object
            longest_contig_num = len(contig.seq) #adding number length of contig to the longest

    #trying to debug
    if longest_contig: #add if there is a valid contig
        long_contig[srr] = longest_contig #putting the longest contig into the dictionary
    else:
        print(f"No valid contigs")

#print(long_contig)


# b) download the betaherpes genomes
betaherpes_download = os.path.join(output_dir, "betaherpes_download") #place to put the final download
os.makedirs(betaherpes_download, exist_ok=True) #mkaing directory if not already made

#downloading betaherpes genomes
subprocess.run(["datasets", "download", "virus", "genome", "taxon", "betaherpesvirinae", "--refseq", "--include", "genome"], check=True)

#find the download so that we can unzip it
downloaded_beta = glob.glob("ncbi_dataset*.zip") #looking for the zip file
if downloaded_beta:
    downloaded_beta = downloaded_beta[0] #getting the first zip file
    unzip_path = os.path.join(betaherpes_download, "betaherpes_download.zip") #making the path which would move the file to a better place

    #moving the file
    os.rename(downloaded_beta, unzip_path)

    #now we can finally unzip it
    with zipfile.ZipFile(unzip_path, "r") as z: #using program to unzip 
        z.extractall(betaherpes_download) #extracting everything in this pathway
    print(f"extraction worked YAY")
else:
    print(f"extraction failed :( )")


# c) make the local database

genomic_fna = os.path.join(output_dir, "betaherpes_download/ncbi_dataset/data/genomic.fna") #place of fasta file
betaherpes_database = os.path.join(output_dir, "Betaherpes_database/betaherpes_database") #output place to put the database

#running the command to make the database
subprocess.run(["makeblastdb", "-in", genomic_fna, "-out", betaherpes_database, "-title", "Betaherpes", "-dbtype", "nucl"], check = True)


'''
#ignore getting rid of
full_conditions = "/home/2025/aavalos4/PythonPipeline_Lexi_Avalos/full_conditions.txt"
log_file = "/home/2025/aavalos4/PythonPipeline_Lexi_Avalos/PipelineProject.log"
'''


#getting just the donor part to add into the final log output
donor_extract = {} #empty dictionary to store the srr associated with each donor

#get donor info
with open(full_conditions, "r") as f: #reading the file
    for line in f: #for every line in the file
        sections = line.strip().split(": ") #split at the colon and take out whitespace 
        if len(sections) == 2: #it should cuz only 2 parts
            donor, srr = sections #get donor and srr number
            donor = " ".join(donor.split()[:2]) #getting just the Donor 1 or Donor 2
            donor_extract[srr] = donor # adding to dictionary as srr as key and donor as value

print(donor_extract)

#getting the longest contigs for the donor rather than each SRR
donor_long_contig = {} #empty dictionary to put the longest contig

for srr, contig in long_contig.items(): #go over the longest contigs for each SRR
    donor_use = donor_extract.get(srr, "Unknown donor") #get the owner or put default as none
    if donor_use not in donor_long_contig: #if the donor is not already in the long contig dictionary
        donor_long_contig[donor_use] = [] # make list for this donor

    donor_long_contig[donor_use].append(contig) #add the longest contig

print(donor_long_contig)

# d) now running blast for each longest contig
with open(log_file, "a") as output: # to add to the log file
    for donor, contig in donor_long_contig.items(): #for each SRR and donor longest contigs in the dictionary
        temp_file = os.path.join(output_dir, f"{donor}_temp.fasta") #making a temporary file to use
        with open(temp_file, "w") as f: #opening the file to write into it 
            SeqIO.write(contig, f, "fasta") #writing to the file

        #making an output for the blast to go
        blast_output = os.path.join(output_dir, f"{donor}_blast_output.tsv")
        
        #hold command for debugging
        blast_command = ["blastn", "-query", temp_file, "-db", betaherpes_database, 
                        "-out", blast_output, "-outfmt", "6 sacc pident length qstart qend sstart send bitscore evalue stitle"] #the tab option with the categories we want
                
        
        #running blast
        subprocess.run(blast_command, check=True) 

        #debug 
        #print(f"running Blast: {' '.join(blast_command)}")

        #write the blast results to the log file
    
        output.write(f"{donor}:\n") #writing the donor name and then next will be the categories
        output.write("sacc\tpident\tlength\tqstart\tqend\tsstart\tsend\tbitscore\tevalue\tstitle\n") #writing all of the categories and separating by a tab
        #print(blast_output)
        #only want the top 10 hits
        with open(blast_output, "r") as b_result:
            for i in b_result.readlines()[:10]: #reading the first 10 lines
                output.write(i) #writing the line that is in the top 10
            #print(i)
        output.write("\n") #putting a space for the next donor

    







            




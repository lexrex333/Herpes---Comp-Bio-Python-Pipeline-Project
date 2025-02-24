# Herpes---Comp-Bio-Python-Pipeline-Project
Hi, welcome to a walkthrough through our 2025 Computational Biology Python Pipeline Project.

This project revolves around the human Herpesvirus 5, which is also known as the human cytomegalovirus (HCMV). 

Throughout this project, we will be using the HCMV transcriptomes of 2 patient donors: 2 transcriptomes each (one consisting from 2 days post infection (dpi), and the other from 6 days post infection) - so a total of 4 transcriptomes! 

If you are curious about the study that sequenced these transcriptomes: [Study](https://pubmed.ncbi.nlm.nih.gov/29158406/).

I have added [sample data](https://github.com/lexrex333/Herpes---Comp-Bio-Python-Pipeline-Project/tree/main/fastq) that is already in fastq format and in a smaller size (10,000 reads) that can be used instead of doing the steps below to get the fastqs. This way it is quicker and more efficient to see if the pipeline works :) ! Please look at section (Step 2: Using sample data) to get help on how to use it. 
## Step 1: Get the Data!
#### 1) Go to the SRA page for each transcriptome:

>>>> a. [Donor 1 - 2dpi](https://www.ncbi.nlm.nih.gov/sra/SRX2896360) 

>>>> b. [Donor 1 - 6dpi](https://www.ncbi.nlm.nih.gov/sra/SRX2896363) 

>>>> c. [Donor 3 - 2dpi](https://www.ncbi.nlm.nih.gov/sra/SRX2896374)

>>>> d. [Donor 3 - 6dpi](https://www.ncbi.nlm.nih.gov/sra/SRX2896375)

You will come across a page that looks like this: 
![image](https://github.com/user-attachments/assets/b9904656-ab7e-4a88-bdde-32a3814db863)

You want to click on the SRR number located under Runs:
![image](https://github.com/user-attachments/assets/fdacf67c-9c99-4a41-b1b9-18b6a8e6d386)

It will then take you to a page that looks like this:
![image](https://github.com/user-attachments/assets/5e15a5d5-af6b-430b-ad58-f92d2ebaa85b)

Click on the Data Access tab:

![image](https://github.com/user-attachments/assets/c20a43a6-6556-4e06-8957-b632c704a413)

Copy the first link with the SRA Normalized Type and AWS Location:
![image](https://github.com/user-attachments/assets/543df29d-2648-46db-a66c-f695591e7897)

#### 2) Download the SRR files
We will be using VSCode and our remote server to download the files through the command line.

So open up the remote server:
![image](https://github.com/user-attachments/assets/d8c98f75-b92a-4b7e-be09-38dc1e42243f)

Optional: Make a directory to put all of your information to be more organized:
```bash
mkdir PipelineProject_Lexi_Avalos 
```
(If you decide to do this, make sure to cd into the directory when downloading the files.)


Go to the terminal and type the command - keep in mind that the URL is the one you copied in the previous step, which will be a new URL for each transcriptome: 
```bash
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660030/SRR5660030
```

When downloading, it should look something like this:
![image](https://github.com/user-attachments/assets/148f23d9-2722-4885-a4bd-111134779f14)

#### 3) Repeat steps 1 and 2 for the rest of the transcriptomes
Here are the rest of the commands for easier access:
```bash
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660033/SRR5660033
```
```bash
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660044/SRR5660044
```
```bash
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660045/SRR5660045
```
Once all are done downloading, your explorer should contain all of them:
![image](https://github.com/user-attachments/assets/499f1074-446f-40bc-8e74-7afe32e8ef68)

#### 4) Changing the SRR file into FASTQ format using Fasterq-dump
Now we must change the SRR files we have, into FASTQ format to uncompress all the data we have.

Make sure you have [fasterq-dump](https://rnnh.github.io/bioinfo-notebook/docs/fasterq-dump.html) installed. 

Now, we are going to use a command for each SRR file (make sure you are in the same directory as before):
```bash
fasterq-dump SRR5660030
```
```bash
fasterq-dump SRR5660033
```
```bash
fasterq-dump SRR5660044
```
```bash
fasterq-dump SRR5660045
```

It will output 2 fastq files for each SRR file:

![image](https://github.com/user-attachments/assets/856f1b9e-fed8-4bba-af01-68a7d4d02658)

After doing that to all SRR files, you now have the fastq files you will use for the rest of this project!!

Optional: Organize the fastq files and SRR files into their own directory within the project directory to be more organized.

![image](https://github.com/user-attachments/assets/01f76de7-a6eb-4352-a985-179bcc501dfe)


## Step 2: Using Sample Data
You will need to make sure to do this to run this pipeline:
#### 1. Make a directory to put all of your data:
For example, for me, it is PythonPipeline_Lexi_Avalos:![image](https://github.com/user-attachments/assets/04cf7abb-cff0-4470-aa69-5fea8c153958). 

This is necessary because you will need to put all of your downloaded data in this specific directory.

#### 2. Make sure everything you need is downloaded:
1. [fastq](https://github.com/lexrex333/Herpes---Comp-Bio-Python-Pipeline-Project/tree/main/fastq)
   
2. [SRR](https://github.com/lexrex333/Herpes---Comp-Bio-Python-Pipeline-Project/tree/main/SRR_file) (These are empty, but you need them to run the samples. If you were running your own data, you would already have these SRR number files in your directory.)

3. [conditions.txt](https://github.com/lexrex333/Herpes---Comp-Bio-Python-Pipeline-Project/blob/main/conditions.txt)
  
4. [full_conditions.txt](https://github.com/lexrex333/Herpes---Comp-Bio-Python-Pipeline-Project/blob/main/full_conditions.txt)
 
5. [sleuth.R](https://github.com/lexrex333/Herpes---Comp-Bio-Python-Pipeline-Project/blob/main/sleuth.R)

6. [Python Wrapper](
   
#### 3. Add all downloaded data into your directory:
I personally did this manually as I dragged the files from my file explorer and put them in my VSCode explorer, directly into the PythonPipeline_Lexi_Avalos directory. It should end up looking like this:

![image](https://github.com/user-attachments/assets/8fc47d18-6369-44e0-8813-9253760512f6)

#### 4. Change the output directory at the top of each wrapper:
Within the [Python Wrapper] : Change the output_dir to match yours:

![image](https://github.com/user-attachments/assets/0bea18bb-c278-4b54-a0f1-d97eee42d0dc)

Within the [Sleuth R Script](https://github.com/lexrex333/Herpes---Comp-Bio-Python-Pipeline-Project/blob/main/sleuth.R): Change the output_dir to match yours:

![image](https://github.com/user-attachments/assets/3274a4b4-fa76-492a-b1d8-bfb4c88dd31c)

#### 5. Dependencies
You will need to make sure you have the following working to use this wrapper, a lot of these come with Python when you download [Python](https://www.python.org/downloads/): [sys](https://www.geeksforgeeks.org/python-sys-module/), [os](https://docs.python.org/3/library/os.html), [subprocess](https://www.geeksforgeeks.org/python-subprocess-module/), [numpy](https://numpy.org/), [pandas](https://pandas.pydata.org/), [zipfile](https://www.geeksforgeeks.org/working-zip-files-python/), [glob](https://docs.python.org/3/library/glob.html), [shutil](https://docs.python.org/3/library/shutil.html), [SeqIO](https://biopython.org/wiki/SeqIO), and [Entrez](https://biopython.org/docs/1.76/api/Bio.Entrez.html). 

#### 6. RUN IT!
Using only one command line, it should run:
```bash
python3 /home/2025/aavalos4/Herpes---Comp-Bio-Python-Pipeline-Project/Pipeline_analysis.py
```
P.S.: You will need to change the pathway to the pathway that you put the python script in to use it. 

#### 7. Examine your log file
YAY - Hopefully it worked, and you should have a log file that looks like this: 







 








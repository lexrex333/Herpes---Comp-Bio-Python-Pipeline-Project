# Herpes---Comp-Bio-Python-Pipeline-Project
Hi, welcome to a walkthrough through our 2025 Computational Biology Python Pipeline Project.

This project revolves around the human Herpesvirus 5, which is also known as the human cytomegalovirus (HCMV). 

Throughout this project, we will be using the HCMV transcriptomes of 2 patient donors: 2 transcriptomes each (one consisting from 2 days post infection (dpi), and the other from 6 days post infection) - so a total of 4 transcriptomes! 

If you are curious about the study that sequenced these transcriptomes: [Study](https://pubmed.ncbi.nlm.nih.gov/29158406/).

## Step 1: Get Raw Data
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

#### 2) Download the .fastq files
We will be using VSCode and our remote server to download the files through the command line.

So open up the remote server:
![image](https://github.com/user-attachments/assets/d8c98f75-b92a-4b7e-be09-38dc1e42243f)

Optional: Make a directory to put all of your information to be more organized:
```bash
mkdir Python_Pipeline_Project
```

Go to the terminal and type the command - keep in mind that the URL is the one you copied in the previous step, which will be a new URL for each transcriptome: 
```bash
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660030/SRR5660030
```

#### 3) Repeat steps 1 and 2 for the rest of the transcriptomes


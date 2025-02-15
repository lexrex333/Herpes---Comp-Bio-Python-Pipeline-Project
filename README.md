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

#### 2) 

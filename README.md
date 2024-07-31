# Opt.-selection-of-negatives
 Optimizing the selection of negative examples for computational RBP binding site prediction

## Table of Contents

- [Project Title](#project-title)
- [Description](#description)
- [Installation](#installation)
- [Usage](#usage)
- [Contributing](#contributing)
- [License](#license)
- [Contact](#contact)

## Description

This project focuses on the impact of the selection of the negative examples to fine-tune a transformer model.
Predicting RNA-binding protein (RBP) binding sites using the DNABERT transformer, specifically tailored for 3-mers. The model is fine-tuned for a binary classification task, distinguishing whether a DNA sequence of variable length is a binding site or not.
We tried three distinct approaches to fine-tune the DNABERT model:

1. **Equal Positive-Negative Ratio (1:1):**
   - In this approach, we used a balanced ratio of 1:1 for positive (binding site) to negative (non-binding site) sequences during the fine-tuning process.

2. **Imbalanced Ratio (1:3):**
   - Here, we experimented with an imbalanced ratio of 1:3, introducing three times as many negative sequences compared to positive ones during the fine-tuning stage.

3. **Selective Similar Negatives (1:1):**
   - This approach involved maintaining a 1:1 ratio but selecting negative sequences that were particularly similar to positive sequences. The goal was to observe the model's performance when dealing with closely related negatives.

## Installation

1. **Install Miniconda:**
   - Download and install Miniconda3 from the [official website](https://docs.anaconda.com/free/miniconda/).

2. **Setup Conda Environment:**
   ```bash
   # Clone the repository
   git clone https://github.com/Luke-hey/Opt.-selection-of-negatives.git
   
   # Navigate to the project directory
   cd Opt.-selection-of-negatives

   # Create and activate Conda environment
   conda create --name your_env_name python=3.8 --file requirements.txt
   conda activate your_env_name


## Usage

### 1. Reference Genome and Annotation

- **Reference Genome (hg38):**
  - Download the reference genome from [UCSC Genome Browser](https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/)
  - Move the reference genome into the root directory of the repository. (/Opt.-selection-of-negatives)
  - unzip the file with gunzip (optional)

- **Annotation of Genomic Regions:**
  - Download the genomic region annotation file from [GENCODE](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz)
  - Same as the reference genome file

### 2. Set up the directory structure

- **Run setup script:**
  - Use following script in the root directory:
    ```bash
    ./setup_scripts/setup_protein_dirs.sh 3 protein1 protein2
    ```
  - This will create a new directory "proteins" where all the single proteins directorys will be saved. Where the 3 is the number of runs, and protein 1 and 2 are the name of the proteins.
  - At the moment the pipeline is working for 3 Run directorys.
    
### 3. Load Positives

- **Download Positives:**
  - Obtain high confidence and low confidence positive sequences from [ENCODE Project](https://www.encodeproject.org/)
  - Note: we are using eclip in this project and we for the low confidence and high confidence positives we are using the bednarrow peaks with the 1 and 2 and 1,2 replicates.
  - This step is important for the name convention of the positives: 1 and 2 are the low confidence positives and 1,2 are the high confidence positives.
  - example names for the 3 files: proteinname_celltype_replicate.bed so for example QKI_K562_1.bed, QKI_K562_2.bed, QKI_K562_1_2.bed
  - move them in the respective protein directory

### 4. Sample Negatives

- **Run sample script:**
  - Use following script in the respective protein directory (/Opt.-selection-of-negatives/proteins/proteinname) to sample negatives for that protein:
    ```bash
    ../../setup_scripts/run_sampling.sh -p proteinname -c celltype -n number of negatives (multiplyer of positives)
    ```
  - TODO: implement queue system to sample negatives for multiple proteins
  - Note: the procedure to obtain the negatives we did the following: 1. sample negatives with n=3, 2. pick 1/3 of the 3x negative set randomly to obtain the 1x negative set, 3. to get the similar negatives we used the tool "fasta-shuffle-letters" with:
   ```bash
    fasta-shuffle-letters -k 2 -dna proteinname_train.fa
    ```
 
### 5. Finetune DNABERT

- **Run DNABERT Script:**
  - Use following script in the proteins directory (/Opt.-selection-of-negatives/proteins) to start the finetuning process:
    ```bash
    ../setup_scripts/run_DNABERT.sh proteinname1 proteinname2
    ```
  - Note: this will start a queue

 ### 6. Prediction on gene regions

- **Run pred_genomic_region Script:**
  - Important: gene regions and test positives should be in bed format, additionally they should be of the format proteinname_sites.bed for positives and proteinname_regions.bed and should be in the respective    protein directory (/Opt.-selection-of-negatives/proteins/proteinname)
  - Use following script in the proteins directory (/Opt.-selection-of-negatives/proteins) to start the predictions:
    ```bash
    ../setup_scripts/run_plot.sh proteinname1 proteinname2
    ```
  - This again will start again a queue
  - Additionally it will produce a plot with the average AUPRC and variance of all the 3 Run directorys

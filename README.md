# Opt.-selection-of-negatives
 Optimizing the selection of negative examples for computational RBP binding site prediction

## Table of Contents

- [Project Title](#project-title)
- [Description](#description)
- [Features](#features)
- [Installation](#installation)
- [Usage](#usage)
- [Contributing](#contributing)
- [License](#license)
- [Contact](#contact)

## Description

This project focuses on the impact of the selection of the negative examples to fine-tune a transformer model.
Predicting RNA-binding protein (RBP) binding sites using the DNABERT transformer, specifically tailored for 3-mers. The model is fine-tuned for a binary classification task, distinguishing whether a DNA sequence of length 101 is a binding site or not.
We tried three distinct approaches to fine-tune the DNABERT model:

1. **Equal Positive-Negative Ratio (1:1):**
   - In this approach, we used a balanced ratio of 1:1 for positive (binding site) to negative (non-binding site) sequences during the fine-tuning process.

2. **Imbalanced Ratio (1:3):**
   - Here, we experimented with an imbalanced ratio of 1:3, introducing three times as many negative sequences compared to positive ones during the fine-tuning stage.

3. **Selective Similar Negatives (1:1):**
   - This approach involved maintaining a 1:1 ratio but selecting negative sequences that were particularly similar to positive sequences. The goal was to observe the model's performance when dealing with closely related negatives.

## Features

key features

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

### 2. Load Positives

- **Download Positives:**
  - Obtain high confidence and low confidence positive sequences from [ENCODE Project](https://www.encodeproject.org/)
  - Note: we are using eclip in this project and we for the low confidence and high confidence positives we are using the bednarrow peaks with the 1 and 2 and 1,2 replicates.
  - This step is important for the name convention of the positives: 1 and 2 are the low confidence positives and 1,2 are the high confidence positives.
  - example names for the 3 files: proteinname_celltype_replicate.bed so for example QKI_K562_1.bed, QKI_K562_2.bed, QKI_K562_1_2.bed

### 3. Sample Negatives

- **Run Main Script:**
  - Use `main_script.py` to sample negatives based on high confidence and both low confidence BED files.
  - Use following command to check input parameters for script:
    ```bash
    python main_script.py --help
    ```
### 4. Finetune DNABERT

- **Run DNABERT_3mer Script:**
  - Use `DNABERT_3mer.py` to finetune the model for the binary classification task using the high confident positives and sampled negatives in fasta format.
  - Use following command to check input parameters for script:
    ```bash
    python DNABERT_3mer.py --help
    ```

### 5. Prediction on test set

- **Run pred_test_seq Script:**
  - After Step 4 there will be a file "test_positive_sequences.fa" which contains the sequences from the test set after the finetuning.
  - Now need to sample some new negatives which weren't in the training set (use generate_negatives.py or sampling_hardneg.py manually)
  - After that use `pred_test_seq.py` to calculate F1 score on the test positive and new sampled negatives which don't overlap with training negatives.
  - Use following command to check input parameters for script:
    ```bash
    python pred_test_seq.py --help
    ```

 ### 6. Prediction on gene regions

- **Run pred_genomic_region Script:**
  - Here we need the genomic region and the positive sequences on those regions in BED format which should be predicted.
  - Might need to pad the positives before using them with the pad_positive.py script
  - Additionally need to finetune DNABERT again but this time check for overlapp between high_conf_positives and your predicted sites. (use bedtools substract to modify the positive file)
  - Negatives might need to be adjusted in size after the positives got modified.
  - After that use `pred_genomic_region.py` to calculate AUPRC and plot the results.
  - Important note: the order of the finetuned model as input of the script should be 1:1 / 1:3 / 1:1 similar / 1:3 mixed
  - Use following command to check input parameters for script:
    ```bash
    python pred_genomic_region.py --help
    ```

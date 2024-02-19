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

1. Install miniconda3 from the official website: link
2. 

```bash
# Example installation commands
git clone https://github.com/yourusername/your-repository.git
cd your-repository
npm install
````

## Usage

import argparse
from Bio import SeqIO
import numpy as np
import pandas as pd
import random
import torch
from datasets import Dataset, load_metric
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import train_test_split
#from torch.utils.data import Dataset
import os
from transformers import (
    AutoTokenizer,
    AutoModelForSequenceClassification,
    AutoConfig,
    DataCollatorWithPadding,
    TrainingArguments,
    Trainer,
    EarlyStoppingCallback,
)

def load_fasta_sequences(file_path):
    sequences = {}
    with open(file_path, 'r') as fasta_file:
        for record in SeqIO.parse(fasta_file, 'fasta'):
            header = record.id
            sequence = str(record.seq).replace('U', 'T')
            sequences[header] = sequence.upper()
    return sequences

def main():
    parser = argparse.ArgumentParser(description="Train DNA BERT for binary classification.")
    parser.add_argument("--positive_sequences", required=True, help="Path to the positive sequences FASTA file")
    parser.add_argument("--negative_sequences", required=True, help="Path to the negative sequences FASTA file")
    args = parser.parse_args()


    negative_sequences_file_name = os.path.splitext(os.path.basename(args.negative_sequences))[0]

    sequences_dict_pos = load_fasta_sequences(args.positive_sequences)
    sequences_dict_neg = load_fasta_sequences(args.negative_sequences)


    # initialize model and tokenizer
    tokenizer = AutoTokenizer.from_pretrained("zhihan1996/DNA_bert_3", use_fast=True, trust_remote_code=True)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print("Using device:", device)
    model = AutoModelForSequenceClassification.from_pretrained("zhihan1996/DNA_bert_3", trust_remote_code=True, num_labels=2).to(device)
    data_collator = DataCollatorWithPadding(tokenizer=tokenizer)


    # needed for tokenization
    def kmer(seq, K = 3):
        kmer_list = []
        kmer_string = ""
        for x in range(len(seq) - K + 1):
            kmer_list.append(seq[x:x+K].upper())
        for i in range(len(kmer_list)):
          if i < len(kmer_list) - 1:
            kmer_string += kmer_list[i] + " "
          else:
            kmer_string += kmer_list[i]
        return kmer_string


    # creating datasets (training 70% train and val are 15% each)
    sequences_pos = list(sequences_dict_pos.values())
    sequences_neg = list(sequences_dict_neg.values())

    # Apply kmer to all sequences
    sequences_pos = [kmer(seq) for seq in sequences_pos]
    sequences_neg = [kmer(seq) for seq in sequences_neg]

    # Concatenate positive and negative sequences
    all_sequences = sequences_pos + sequences_neg

    # Labels for positive and negative sequences
    labels = [1] * len(sequences_pos) + [0] * len(sequences_neg)

    # Set the random seed for reproducibility
    seed = 42

    # Print the class distribution before splitting
    unique_classes, counts = np.unique(labels, return_counts=True)
    print("Class distribution before splitting:", dict(zip(unique_classes, counts)))

    # Split into train (70%), validation (15%), and test (15%) sets
    train_sequences, test_val_sequences, train_labels, test_val_labels = train_test_split(
        all_sequences, labels, test_size=0.3, random_state=seed, stratify=labels
    )
    # Print the class distribution after the first split
    unique_classes_train, counts_train = np.unique(train_labels, return_counts=True)
    print("Class distribution in training set:", dict(zip(unique_classes_train, counts_train)))

    # Split test_val_sequences into validation (50%) and test (50%) sets
    val_sequences, test_sequences, val_labels, test_labels = train_test_split(
        test_val_sequences, test_val_labels, test_size=0.5, random_state=seed, stratify=test_val_labels
    )
    # Print the class distribution after the second split
    unique_classes_val, counts_val = np.unique(val_labels, return_counts=True)
    print("Class distribution in validation set:", dict(zip(unique_classes_val, counts_val)))
    # Create dictionaries for each split
    def create_dict(sequences, labels, tokenizer):
        encodings = tokenizer(sequences, padding=True, return_tensors="pt")
        return {
            'input_ids': encodings['input_ids'].tolist(),
            'token_type_ids': encodings['token_type_ids'].tolist(),
            'attention_mask': encodings['attention_mask'].tolist(),
            'labels': labels
        }

    train_dict = create_dict(train_sequences, train_labels, tokenizer)
    val_dict = create_dict(val_sequences, val_labels, tokenizer)
    test_dict = create_dict(test_sequences, test_labels, tokenizer)

    # Create datasets
    train_ds = Dataset.from_dict(train_dict)
    val_ds = Dataset.from_dict(val_dict)
    test_ds = Dataset.from_dict(test_dict)

    if not os.path.exists("test_positive_sequences.fa"):
        test_positive_sequences = [tokenizer.decode(test_ds['input_ids'][i], skip_special_tokens=True) for i, label in enumerate(test_ds['labels']) if label == 1]

        # needed to save positives in test set for predictions
        def save_sequences_to_fasta(sequences, file_path):
            with open(file_path, 'w') as fasta_file:
                for i, seq in enumerate(sequences):
                    fasta_file.write(f'>positive_sequence_{i}\n{seq}\n')

        # Save the positive sequences to a file / need to do this only once per dataset
        positive_sequences_file_path = "test_positive_sequences.fa"
        save_sequences_to_fasta(test_positive_sequences, positive_sequences_file_path)




    def compute_metrics(p):
        predictions = p.predictions
        labels = p.label_ids
        preds = torch.softmax(torch.tensor(predictions), dim=1)[:, 1]
        roc_auc = roc_auc_score(labels, preds)
        return {"roc_auc": roc_auc}


    batch_size = 16
    training_args = TrainingArguments(output_dir=f"DNABERT3mer_results_{negative_sequences_file_name}",
                                      evaluation_strategy = "epoch",
                                      #eval_steps = 100,
                                      save_strategy = "epoch",
                                      num_train_epochs=100,
                                      learning_rate=2e-5,
                                      per_device_train_batch_size=batch_size,
                                      per_device_eval_batch_size=batch_size,
                                      weight_decay=0.01,
                                      metric_for_best_model = "eval_loss",
                                      load_best_model_at_end=True,
                                      #push_to_hub=True
                                      )

    trainer = Trainer(
        model,
        training_args,
        train_dataset=train_ds,
        eval_dataset=val_ds,
        data_collator=data_collator,
        tokenizer=tokenizer,
        compute_metrics=compute_metrics,
        callbacks = [EarlyStoppingCallback(early_stopping_patience=6)]
    )

    trainer.train()
    results = trainer.evaluate(test_ds)

    print(results)


    # Save the model and tokenizer with a name based on the negative sequences file
    save_model_path = f"finetuned_DNABERT3mer_{negative_sequences_file_name}"
    save_tokenizer_path = f"finetuned_DNABERT3mer_{negative_sequences_file_name}"

    # Uncomment the following lines to save the model and tokenizer
    trainer.save_model(save_model_path)
    tokenizer.save_pretrained(save_tokenizer_path)


if __name__ == "__main__":
    main()

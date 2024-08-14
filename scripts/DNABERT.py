import argparse
from Bio import SeqIO
import numpy as np
import pandas as pd
import random
import torch
from datasets import Dataset, load_metric
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import train_test_split
import shutil
import os
from transformers import (
    AutoTokenizer,
    AutoModelForSequenceClassification,
    DataCollatorWithPadding,
    TrainingArguments,
    Trainer,
    EarlyStoppingCallback,
)
import torch.optim as optim
os.environ["TOKENIZERS_PARALLELISM"] = "false"

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
    parser.add_argument("--output_dir", default=".", help="Output directory to save the model and tokenizer")
    parser.add_argument("--kmer_length", type=int, choices=[3, 4, 5, 6], default=3, help="Length of k-mer to use (default: 3)")
    args = parser.parse_args()

    # Set the random seed for reproducibility
    seed = 42
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)

    negative_sequences_file_name = os.path.splitext(os.path.basename(args.negative_sequences))[0]

    sequences_dict_pos = load_fasta_sequences(args.positive_sequences)
    sequences_dict_neg = load_fasta_sequences(args.negative_sequences)

    # initialize model and tokenizer
    model_name = f"zhihan1996/DNA_bert_{args.kmer_length}"
    tokenizer = AutoTokenizer.from_pretrained(model_name, use_fast=True, trust_remote_code=True)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print("Using device:", device)
    model = AutoModelForSequenceClassification.from_pretrained(model_name, trust_remote_code=True, num_labels=2).to(device)
    data_collator = DataCollatorWithPadding(tokenizer=tokenizer)

    # Define kmer function
    def kmer(seq, K=3):
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

    sequences_pos = list(sequences_dict_pos.values())
    sequences_neg = list(sequences_dict_neg.values())

    # Apply kmer to all sequences
    sequences_pos = [kmer(seq, args.kmer_length) for seq in sequences_pos]
    sequences_neg = [kmer(seq, args.kmer_length) for seq in sequences_neg]

    # Concatenate positive and negative sequences
    all_sequences = sequences_pos + sequences_neg

    # Labels for positive and negative sequences
    labels = [1] * len(sequences_pos) + [0] * len(sequences_neg)

    # Print the class distribution before splitting
    unique_classes, counts = np.unique(labels, return_counts=True)
    print("Class distribution before splitting:", dict(zip(unique_classes, counts)))

    # Split into train (80%) and validation (20%) sets
    train_sequences, val_sequences, train_labels, val_labels = train_test_split(
        all_sequences, labels, test_size=0.2, random_state=seed, stratify=labels
    )

    # Print the class distribution in training and validation sets
    unique_classes_train, counts_train = np.unique(train_labels, return_counts=True)
    print("Class distribution in training set:", dict(zip(unique_classes_train, counts_train)))

    unique_classes_val, counts_val = np.unique(val_labels, return_counts=True)
    print("Class distribution in validation set:", dict(zip(unique_classes_val, counts_val)))

    # Create dictionaries for each split
    def create_dict(sequences, labels, tokenizer):
        encodings = tokenizer(sequences, padding=True, truncation=False, return_tensors="pt")
        return {
            'input_ids': encodings['input_ids'].tolist(),
            'token_type_ids': encodings['token_type_ids'].tolist(),
            'attention_mask': encodings['attention_mask'].tolist(),
            'labels': labels
        }

    train_dict = create_dict(train_sequences, train_labels, tokenizer)
    val_dict = create_dict(val_sequences, val_labels, tokenizer)

    # Create datasets
    train_ds = Dataset.from_dict(train_dict)
    val_ds = Dataset.from_dict(val_dict)

    def compute_metrics(p):
        predictions = p.predictions
        labels = p.label_ids
        preds = torch.softmax(torch.tensor(predictions), dim=1)[:, 1]
        roc_auc = roc_auc_score(labels, preds)
        return {"roc_auc": roc_auc}

    batch_size = 16
    output_dir = f"DNABERT{args.kmer_length}_results_{negative_sequences_file_name}"
    training_args = TrainingArguments(
        output_dir=output_dir,
        evaluation_strategy="epoch",
        save_strategy="epoch",
        num_train_epochs=100,
        learning_rate=2e-5,
        per_device_train_batch_size=batch_size,
        per_device_eval_batch_size=batch_size,
        weight_decay=0.01,
        metric_for_best_model="eval_loss",
        load_best_model_at_end=True,
        seed=seed,
    )

    trainer = Trainer(
        model,
        training_args,
        train_dataset=train_ds,
        eval_dataset=val_ds,
        data_collator=data_collator,
        tokenizer=tokenizer,
        compute_metrics=compute_metrics,
        callbacks=[EarlyStoppingCallback(early_stopping_patience=6)],
        optimizers=(optim.AdamW(model.parameters(), lr=2e-5), None)
    )

    trainer.train()

    # Save the model and tokenizer
    save_model_path = os.path.join(args.output_dir, f"finetuned_DNABERT{args.kmer_length}_{negative_sequences_file_name}")
    save_tokenizer_path = os.path.join(args.output_dir, f"finetuned_DNABERT{args.kmer_length}_{negative_sequences_file_name}")
    trainer.save_model(save_model_path)
    tokenizer.save_pretrained(save_tokenizer_path)

    # Clean up output directory after training (if necessary)
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
        print(f"Removed directory: {output_dir}")

if __name__ == "__main__":
    main()

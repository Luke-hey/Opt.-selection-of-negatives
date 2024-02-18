import pybedtools
import torch
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math
import argparse

from Bio import SeqIO
from Bio.Seq import Seq

def load_fasta_sequences(file_path):
    sequences = {}
    with open(file_path, 'r') as fasta_file:
        for record in SeqIO.parse(fasta_file, 'fasta'):
            header = record.id
            sequence = str(record.seq).replace('U', 'T')
            if len(sequence) >= 101:
                sequences[header] = sequence.upper()
    return sequences

def normalize_vector(vec):
    total = sum(vec)
    if total == 0:
        return vec
    return [x / total for x in vec]

def dinucleotide_vector(dna_sequence):
    dinucleotides = ['AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT']
    counts = [dna_sequence.count(dinuc) for dinuc in dinucleotides]
    return normalize_vector(counts)

def dict_to_dinucleotide_vectors(dna_dict):
    dinucleotide_dict = {}
    for header, sequence in dna_dict.items():
        dinucleotide_dict[header] = dinucleotide_vector(sequence)
    return dinucleotide_dict

def cosine_similarity(vector_a, vector_b):
    dot_product = np.dot(vector_a, vector_b)
    norm_a = np.linalg.norm(vector_a)
    norm_b = np.linalg.norm(vector_b)
    if norm_a == 0 or norm_b == 0:
        return np.nan

    return dot_product / (norm_a * norm_b)


def generate_hard_negatives(genes_without_pos_file, positives_file, output_file, num_sequences_to_save=4186):
    genes_without_pos = load_fasta_sequences(genes_without_pos_file)
    sequences_dict_pos = load_fasta_sequences(positives_file)

    dinucleotide_vectors = dict_to_dinucleotide_vectors(sequences_dict_pos)

    average_vector = np.zeros(16)

    num_vectors = len(dinucleotide_vectors)
    for vec in dinucleotide_vectors.values():
        average_vector += np.array(vec)

    average_vector /= num_vectors

    # Create a list to store windows that meet the dynamic threshold
    matching_windows = []

    # Define window size and step size
    window_size = 101
    step_size = 1

    # Sort sequences by length in descending order
    sorted_genes_without_pos = sorted(genes_without_pos.items(), key=lambda x: len(x[1]), reverse=True)

    # Iterate over sequences in genes_without_pos
    for header, sequence in sorted_genes_without_pos:
        # Check if the sequence is long enough for a window
        if len(sequence) < window_size:
            continue

        # Iterate over the sequence with the specified step size
        for i in range(0, len(sequence) - window_size + 1, step_size):
            # Extract the window
            window_sequence = sequence[i:i+window_size]

            # Calculate dinucleotide vector for the window
            window_vector = dinucleotide_vector(window_sequence)

            # Calculate cosine similarity with the average vector
            similarity = cosine_similarity(window_vector, average_vector)

            # Check if similarity is above the dynamic threshold
            if len(matching_windows) < num_sequences_to_save:
                matching_windows.append((header, i, window_sequence, similarity))
            else:
                # Sort matching_windows by similarity in descending order
                matching_windows.sort(key=lambda x: x[3], reverse=True)

                # Replace the lowest similarity window if the current window has a higher similarity
                if similarity > matching_windows[-1][3]:
                    matching_windows[-1] = (header, i, window_sequence, similarity)

    with open(output_file, 'w') as fastafile:
        for i, window_info in enumerate(matching_windows, 1):
            header = f"hard_neg_{i}"
            sequence = window_info[2]

            # Write header and sequence to the FASTA file
            fastafile.write(f">{header}\n{sequence}\n")

        print(f"Saved {min(i, num_sequences_to_save)} sequences to {output_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate hard negative sequences based on dinucleotide composition.")
    parser.add_argument("genes_without_pos_file", help="FASTA file containing gene sequences without positive samples")
    parser.add_argument("positives_file", help="FASTA file containing positive samples")
    parser.add_argument("output_file", help="Output file for hard negative sequences")
    parser.add_argument("--num_sequences_to_save", type=int, default=4186, help="Number of sequences to save (default: 4186)")

    args = parser.parse_args()

    generate_hard_negatives(args.genes_without_pos_file, args.positives_file, args.output_file, args.num_sequences_to_save)

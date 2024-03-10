import argparse
import torch
from transformers import AutoModelForSequenceClassification, AutoTokenizer
import pybedtools
from pyfaidx import Fasta
from sklearn.metrics import confusion_matrix, precision_score, recall_score
from tqdm import tqdm
import os
from Bio import SeqIO
os.environ["TOKENIZERS_PARALLELISM"] = "false"

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

# Define the kmer function
def kmer(seq, K=3):
    kmer_list = []
    kmer_string = ""
    for x in range(len(seq) - K + 1):
        kmer_list.append(seq[x : x + K].upper())
    for i in range(len(kmer_list)):
        if i < len(kmer_list) - 1:
            kmer_string += kmer_list[i] + " "
        else:
            kmer_string += kmer_list[i]
    return kmer_string

# Define the prediction function
def pred(seq, model, tokenizer, threshold=0.80):
    with torch.no_grad():
        model.eval()
        inputs = tokenizer(kmer(seq), return_tensors="pt").to(device)
        outputs = model(**inputs)
        logits = outputs.logits
        probabilities = torch.nn.functional.softmax(logits, dim=-1)
        prediction = (probabilities[:, 1] > threshold).long()
        return prediction.item()

# Function to perform sliding window prediction on genomic regions
def predict_genomic_regions(genomic_bed, positive_bed, model, tokenizer):
    # Load bed files
    genomic_regions = pybedtools.BedTool(genomic_bed)
    positive_sequences = pybedtools.BedTool(positive_bed)

    # Load genome reference for sequence extraction
    genome_reference = "/home/bubu23/hg38.fa"        # need to change path here
    genome = Fasta(genome_reference)

    # Initialize an empty list to store positive predicted windows
    positive_predictions = []
    # Iterate over genomic regions with sliding window
    for region in tqdm(genomic_regions, desc="Processing genomic regions", unit="region"):
        chrom, start, end, strand = region.chrom, int(region.start), int(region.end), region.strand

        if strand == '+':
            sequence = genome[chrom][start:end].seq
        elif strand == '-':
            # If the strand is negative, you might want to reverse complement the sequence
            sequence = genome[chrom][start:end].reverse.complement.seq

        for window_start in tqdm(range(0, len(sequence) - 101, 20), desc="Processing windows", unit="window", leave=False):
            window_end = window_start + 101
            window_sequence = sequence[window_start:window_end]

            # Make predictions using the provided function
            prediction = pred(window_sequence, model, tokenizer)

            # Check if the prediction is positive and add to the list
            if prediction == 1:
                # Adjust format to include a fifth column for the name (set to ".")
                positive_predictions.append((chrom, window_start + start, window_end + start, ".", 0, strand))

    """
    # Save positive predicted windows to a BED file
    with open("positive_predictions.bed", "w") as bedfile:
        for prediction in positive_predictions:
            bedfile.write("\t".join(map(str, prediction)) + "\n")
    """
    return positive_predictions


# Function to evaluate predictions using precision, recall, and confusion matrix
def evaluate_predictions(ground_truth_bed, predicted_bed):
    ground_truth = pybedtools.BedTool(ground_truth_bed)
    predicted = pybedtools.BedTool(predicted_bed)



    # Intersect predicted and ground truth
    intersection = ground_truth.intersect(predicted, u=True, s=True)


    TP = len(intersection)
    FP = len(predicted) - TP
    FN = len(ground_truth) - TP

    # Calculate Precision, Recall, and F1
    precision = TP / (TP + FP) if (TP + FP) > 0 else 0
    recall = TP / (TP + FN) if (TP + FN) > 0 else 0
    f1 = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0

    # Print the results
    print("Precision:", precision)
    print("Recall:", recall)
    print("F1 Score:", f1)
    return precision, recall, f1




def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Genomic Region Prediction")
    parser.add_argument("--genomic_bed", required=True, help="Path to the genomic bed file")
    parser.add_argument("--positive_bed", required=True, help="Path to the positive sequences bed file")
    parser.add_argument("--model_name", required=True, help="Pretrained model name")
    args = parser.parse_args()
    model_name = args.model_name
    # Load pretrained model and tokenizer
    tokenizer = AutoTokenizer.from_pretrained(model_name, trust_remote_code=True)
    model = AutoModelForSequenceClassification.from_pretrained(model_name, trust_remote_code=True, num_labels=2).to(device)

    # Perform sliding window prediction on genomic regions
    predicted_windows = predict_genomic_regions(args.genomic_bed, args.positive_bed, model, tokenizer)


    # Evaluate predictions
    evaluate_predictions(args.positive_bed, predicted_windows)

if __name__ == "__main__":
    main()

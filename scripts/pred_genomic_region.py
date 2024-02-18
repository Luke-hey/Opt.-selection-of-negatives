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
def pred(seq, model, tokenizer, threshold=0.96):
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
    genome_reference = "/home/bubu23/hg38.fa"
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

    # Save positive predicted windows to a BED file
    with open("positive_predictions.bed", "w") as bedfile:
        for prediction in positive_predictions:
            bedfile.write("\t".join(map(str, prediction)) + "\n")

    return positive_predictions

    # Merge overlapping or adjacent regions to eliminate duplicates
    processed_predictions = pybedtools.BedTool("positive_predictions.bed").sort().merge()


# Function to evaluate predictions using precision, recall, and confusion matrix
def evaluate_predictions(ground_truth_bed, predicted_bed):
    ground_truth = pybedtools.BedTool(ground_truth_bed)
    predicted = pybedtools.BedTool(predicted_bed)
    print("Length of ground_truth:", len(ground_truth))
    print("Length of predicted:", len(predicted))


    # Intersect predicted and ground truth
    intersection = ground_truth.intersect(predicted, u=True, s=True)


    precision = len(list(intersection)) / len(predicted)
    recall = len(list(intersection)) / len(ground_truth)
    f1_score = 2 * (precision * recall) / (precision + recall)

    print(f"Precision: {precision}")
    print(f"Recall: {recall}")
    print(f"F1 Score: {f1_score}")


### TEST ###
def extract_sequences_from_bed(bed_file, hg_fasta="/home/bubu23/hg38.fa"):
    sequence_dict = {}

    bed_file = pybedtools.BedTool(bed_file)
    bed_file = bed_file.sequence(fi=hg_fasta, s=True, nameOnly=True)

    for record in SeqIO.parse(str(bed_file.seqfn), "fasta"):
        index_of_parenthesis = record.id.find("(")
        result = record.id[:index_of_parenthesis] if index_of_parenthesis != -1 else record.id
        sequence_dict[result] = str(record.seq)

    return sequence_dict



def get_gene_ranges(bed_file):
    # Function to extract gene ranges, chromosome, and strand from a BED file using pybedtools
    gene_info = {}
    bed = pybedtools.BedTool(bed_file)

    for interval in bed:
        header = interval.name
        chrom = interval.chrom
        start = interval.start
        end = interval.end
        strand = interval.strand
        gene_info[header] = {'chrom': chrom, 'start': start, 'end': end, 'strand': strand}

    return gene_info


def split_sequences_with_gene_ranges(sequences_dict, bed_file, model, tokenizer, window_size, step_size):
    windows = []
    gene_info = get_gene_ranges(bed_file)

    # Use tqdm to create a progress bar
    for header, sequence in tqdm(sequences_dict.items(), desc="Processing sequences", unit="sequence"):
        if header in gene_info:
            chrom = gene_info[header]['chrom']
            start = gene_info[header]['start']
            end = gene_info[header]['end']
            strand = gene_info[header]['strand']

            seq_len = len(sequence)

            for i in range(0, seq_len - window_size + 1, step_size):
                window_sequence = sequence[i : i + window_size]
                if pred(window_sequence, model, tokenizer) == 0:
                    continue
                else:
                    if strand == '+':
                        window_start = start + i
                        window_end = window_start + window_size
                    elif strand == '-':
                        window_end = end - i
                        window_start = window_end - window_size

                    window_header = f"{chrom} {window_start} {window_end} {header} 0 {strand}"
                    windows.append((window_header, window_sequence))

    return windows

def create_bed_file(data, output_bed_file):
    headers = [header for header, _ in data]
    bed_tool = pybedtools.BedTool("\n".join(headers), from_string=True)
    bed_tool.saveas(output_bed_file)
### TEST ###

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
    #predicted_windows = predict_genomic_regions(args.genomic_bed, args.positive_bed, model, tokenizer)

    ### TEST ###
    sequences = extract_sequences_from_bed(args.genomic_bed)
    windows = split_sequences_with_gene_ranges(sequences, args.genomic_bed, model, tokenizer, window_size=101, step_size=20)
    create_bed_file(windows, output_bed_file = "output.bed")
    evaluate_predictions(args.positive_bed, "output.bed")
    ### TEST ###
    # Evaluate predictions
    #evaluate_predictions(args.positive_bed, "positive_predictions.bed")

if __name__ == "__main__":
    main()

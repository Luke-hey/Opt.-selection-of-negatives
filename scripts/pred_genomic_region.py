import argparse
import torch
from transformers import AutoModelForSequenceClassification, AutoTokenizer
import pybedtools
from pyfaidx import Fasta
from tqdm import tqdm
import os
from sklearn.metrics import precision_recall_curve, auc
import matplotlib.pyplot as plt

os.environ["TOKENIZERS_PARALLELISM"] = "false"

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

def get_file_path(file_name):
    current_dir = os.path.dirname(os.path.realpath(__file__))
    return os.path.join(current_dir, file_name)

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

def pred_proba(seq, model, tokenizer):
    with torch.no_grad():
        model.eval()
        inputs = tokenizer(kmer(seq), return_tensors="pt").to(device)
        outputs = model(**inputs)
        logits = outputs.logits
        probabilities = torch.nn.functional.softmax(logits, dim=-1)
        return probabilities[:, 1].item()

def predict_genomic_regions(genomic_bed, positive_bed, model, tokenizer):
    genomic_regions = pybedtools.BedTool(genomic_bed)
    positive_sequences = pybedtools.BedTool(positive_bed)
    genome_reference = "/home/bubu23/hg38.fa"
    genome = Fasta(genome_reference)

    positive_predictions = []

    for region in tqdm(genomic_regions, desc="Processing genomic regions", unit="region"):
        chrom, start, end, strand = region.chrom, int(region.start), int(region.end), region.strand

        if strand == '+':
            sequence = genome[chrom][start:end].seq
        elif strand == '-':
            sequence = genome[chrom][start:end].reverse.complement.seq

        for window_start in tqdm(range(0, len(sequence) - 101, 5), desc="Processing windows", unit="window", leave=False):
            window_end = window_start + 101
            window_sequence = sequence[window_start:window_end]
            probability = pred_proba(window_sequence, model, tokenizer)
            if strand == '+':
                positive_predictions.append((chrom, window_start + start, window_end + start, ".", str(probability), strand))
            elif strand == '-':
                positive_predictions.append((chrom, end-window_end, end-window_start, ".", str(probability), strand))
    return positive_predictions

def calculate_precision_recall_curve(ground_truth_bed, predicted_windows):
    ground_truth = pybedtools.BedTool(ground_truth_bed)
    true_labels = []
    probabilities = []

    for interval in predicted_windows:
        intersection = ground_truth.intersect([interval], u=True, s=True, f=0.1)
        if len(intersection) > 0:
            true_labels.append(1)
        else:
            true_labels.append(0)
        probabilities.append(float(interval[4]))

    precision, recall, _ = precision_recall_curve(true_labels, probabilities)
    auprc = auc(recall, precision)

    return precision, recall, auprc

def plot_precision_recall_curves(models_data):
    legend_mapping = {
        models_data[0][0]: "Equal (1:1)",
        models_data[1][0]: "Imbalanced (1:3)",
        models_data[2][0]: "Similar (1:1)",
        models_data[3][0]: "Mixed (1:3)"
    }

    for i, (model_name, precision, recall, auprc) in enumerate(models_data):
        legend_name = legend_mapping[model_name]
        plt.plot(recall, precision, label=f'{legend_name} = {auprc:.2f}')

    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title('Precision-Recall Curve')
    plt.legend(title='AUPRC', fontsize='small')
    plt.show()

def main():
    parser = argparse.ArgumentParser(description="Genomic Region Prediction")
    parser.add_argument("--genomic_bed", required=True, help="Path to the genomic bed file")
    parser.add_argument("--positive_bed", required=True, help="Path to the positive sequences bed file")
    parser.add_argument("--model_names", nargs='+', required=True, help="List of pretrained model names")
    args = parser.parse_args()

    models_data = []

    for model_name in args.model_names:
        tokenizer = AutoTokenizer.from_pretrained(model_name, trust_remote_code=True)
        model = AutoModelForSequenceClassification.from_pretrained(model_name, trust_remote_code=True, num_labels=2).to(device)
        predicted_windows = predict_genomic_regions(args.genomic_bed, args.positive_bed, model, tokenizer)
        precision, recall, auprc = calculate_precision_recall_curve(args.positive_bed, predicted_windows)
        models_data.append((model_name, precision, recall, auprc))

    plot_precision_recall_curves(models_data)

if __name__ == "__main__":
    main()

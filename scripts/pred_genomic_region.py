import argparse
import torch
from transformers import AutoModelForSequenceClassification, AutoTokenizer
import pybedtools
from pyfaidx import Fasta
from tqdm import tqdm
import os
from sklearn.metrics import precision_recall_curve, auc, f1_score
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
from Bio import SeqIO

os.environ["TOKENIZERS_PARALLELISM"] = "false"

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

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

def pred_proba(seq, model, tokenizer, kmer_length):
    with torch.no_grad():
        model.eval()
        inputs = tokenizer(kmer(seq, K=kmer_length), return_tensors="pt").to(device)
        outputs = model(**inputs)
        logits = outputs.logits
        # could use sigmod instead of softmax here instead (more efficent)
        probabilities = torch.nn.functional.softmax(logits, dim=-1)
        return probabilities[:, 1].item()

def predict_genomic_regions(genomic_bed, positive_bed, model, tokenizer, kmer_length):
    genomic_regions = pybedtools.BedTool(genomic_bed)
    script_dir = os.path.dirname(os.path.abspath(__file__))
    genome_reference = os.path.join(script_dir, "../hg38.fa")
    genome = Fasta(genome_reference)

    positive_predictions = []
    window_length = 101
    step_size = 5

    for region in tqdm(genomic_regions, desc="Processing genomic regions", unit="region"):
        chrom, start, end, strand = region.chrom, int(region.start), int(region.end), region.strand

        if strand == '+':
            sequence = genome[chrom][start:end].seq
        elif strand == '-':
            sequence = genome[chrom][start:end].reverse.complement.seq

        for window_start in tqdm(range(0, len(sequence) - window_length, step_size), desc="Processing windows", unit="window", leave=False):
            window_end = window_start + window_length
            window_sequence = sequence[window_start:window_end]
            probability = pred_proba(window_sequence, model, tokenizer, kmer_length)
            if strand == '+':
                positive_predictions.append((chrom, window_start + start, window_end + start, ".", str(probability), strand))
            elif strand == '-':
                positive_predictions.append((chrom, end - window_end, end - window_start, ".", str(probability), strand))
    return positive_predictions

def calculate_precision_recall_curve(ground_truth_bed, predicted_windows):
    ground_truth = pybedtools.BedTool(ground_truth_bed)
    true_labels = []
    probabilities = []

    for interval in predicted_windows:
        intersection = ground_truth.intersect([interval], u=True, s=True, f=0.5)
        if len(intersection) > 0:
            true_labels.append(1)
        else:
            true_labels.append(0)
        probabilities.append(float(interval[4]))

    precision, recall, _ = precision_recall_curve(true_labels, probabilities)
    auprc = auc(recall, precision)

    return precision, recall, auprc

def calculate_f1_score(ground_truth_bed, predicted_windows, threshold):
    ground_truth = pybedtools.BedTool(ground_truth_bed)
    true_labels = []
    predictions = []

    for interval in predicted_windows:
        intersection = ground_truth.intersect([interval], u=True, s=True, f=0.5)
        if len(intersection) > 0:
            true_labels.append(1)
        else:
            true_labels.append(0)
        predictions.append(1 if float(interval[4]) > threshold else 0)

    return f1_score(true_labels, predictions)

def calculate_threshold(positive_fasta, model, tokenizer, kmer_length):
    # Predict probabilities on the positive sequences to determine threshold
    predictions = []
    with open(positive_fasta, 'r') as fasta_file:
        for record in SeqIO.parse(fasta_file, 'fasta'):
            sequence = str(record.seq).replace('U', 'T').upper()
            probability = pred_proba(sequence, model, tokenizer, kmer_length)
            predictions.append(probability)

    sorted_predictions = sorted(predictions, reverse=True)
    threshold_index = int(len(sorted_predictions) * 0.50)
    threshold = sorted_predictions[threshold_index]
    return threshold

def plot_precision_recall_curves(models_data, output_name, output_dir):
    """Plot precision-recall curves for different models with dynamic legend mapping."""

    plt.figure(figsize=(11.69, 8.27))  # Set the figure size to A4 landscape

    for model_name, precision, recall, auprc in models_data:
        plt.plot(recall, precision, label=f'{model_name} = {auprc:.2f}')

    plt.xlabel('Recall', fontsize=18)  # Set x-axis label font size
    plt.ylabel('Precision', fontsize=18)  # Set y-axis label font size
    plt.title('Precision-Recall Curve', fontsize=20)  # Set title font size

    # Set legend font sizes
    plt.legend(title='AUPRC', fontsize=14, title_fontsize='16')

    # Set tick parameters with increased font size
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)

    plt.grid(True) # Add a grid for better readability

    plt.tight_layout()  # Adjust layout to fit everything nicely

    output_path = os.path.join(output_dir, f"{output_name}_precision_recall_curve.png")
    plt.savefig(output_path, dpi=300)  # Save the plot with a high resolution for printing

    #plt.show()


def interpolate_precision_recall(precision, recall, common_recall):
    f = interp1d(recall, precision, kind='linear', bounds_error=False, fill_value=(precision[0], precision[-1]))
    return f(common_recall)

def plot_average_auprc(models_data, output_name, output_dir):
    """Plot average precision-recall curves for different models with dynamic legend mapping."""

    plt.figure(figsize=(11.69, 8.27))  # A4 landscape size

    # Define a common recall range for interpolation
    common_recall = np.linspace(0, 1, 1000)

    # Dynamic legend mapping based on model names
    legend_mapping = {}
    for data in models_data.values():
        for model_name, _, _, _ in data:
            if "Neg1x" in model_name:
                legend_mapping[model_name] = "Equal (1:1)"
            elif "Neg3x" in model_name:
                legend_mapping[model_name] = "Imbalanced (1:3)"
            else:
                legend_mapping[model_name] = "Shuffled (1:1)"

    for data in models_data.values():
        interpolated_precisions = []
        auprcs = []

        for model_name, precision, recall, auprc in data:
            # Interpolate precision at common recall points
            interp_precision = interpolate_precision_recall(precision, recall, common_recall)
            interpolated_precisions.append(interp_precision)
            auprcs.append(auprc)

        # Calculate average precision and AUPRC
        precision_avg = np.mean(interpolated_precisions, axis=0)
        auprc_avg = np.mean(auprcs)

        plt.plot(common_recall, precision_avg, label=f'{legend_mapping[model_name]} = {auprc_avg:.2f}')

    plt.xlabel('Recall', fontsize=18)
    plt.ylabel('Precision', fontsize=18)
    plt.title('Precision-Recall Curve', fontsize=20)
    plt.legend(title='AUPRC', fontsize=14, title_fontsize='16')
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.grid(True)
    plt.tight_layout()

    output_path = os.path.join(output_dir, f"{output_name}_average_precision_recall_curve.png")
    plt.savefig(output_path, dpi=300)
    #plt.show()


def plot_variance_auprc(models_data, output_name, output_dir):
    """Plot variance in precision-recall curves for different model groups with dynamic legend mapping."""

    plt.figure(figsize=(11.69, 8.27))  # A4 landscape size

    # Define a common recall range for interpolation
    common_recall = np.linspace(0, 1, 1000)

    # Dynamic legend mapping based on model names
    legend_mapping = {
        "1:1": "Equal (1:1)",
        "1:3": "Imbalanced (1:3)",
        "shuffled": "Shuffled (1:1)"
    }

    for group, data in models_data.items():
        interpolated_precisions = []
        auprcs = []

        for model_name, precision, recall, auprc in data:
            # Interpolate precision at common recall points
            interp_precision = interpolate_precision_recall(precision, recall, common_recall)
            interpolated_precisions.append(interp_precision)
            auprcs.append(auprc)

        if len(interpolated_precisions) == 0:
            continue
        # Convert list to numpy array for easier manipulation
        interpolated_precisions = np.array(interpolated_precisions)

        # Calculate the minimum and maximum precision values at each recall point
        precision_min = np.min(interpolated_precisions, axis=0)
        precision_max = np.max(interpolated_precisions, axis=0)
        precision_avg = np.mean(interpolated_precisions, axis=0)
        auprc_avg = np.mean(auprcs)

        # Plot the filled area between the minimum and maximum precision values
        plt.fill_between(common_recall, precision_min, precision_max, alpha=0.2, label=f'{legend_mapping[group]} Variance')

        # Plot the average precision-recall curve
        plt.plot(common_recall, precision_avg, label=f'{legend_mapping[group]} Average = {auprc_avg:.2f}')

    plt.xlabel('Recall', fontsize=18)
    plt.ylabel('Precision', fontsize=18)
    plt.title('Precision-Recall Curve with Grouped Variance', fontsize=20)
    plt.legend(title='AUPRC', fontsize=14, title_fontsize='16')
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.grid(True)
    plt.tight_layout()

    output_path = os.path.join(output_dir, f"{output_name}_grouped_variance_precision_recall_curve.png")
    plt.savefig(output_path, dpi=300)
    #plt.show()

def main():
    parser = argparse.ArgumentParser(description="Genomic Region Prediction")
    parser.add_argument("--genomic_bed", required=True, help="Path to the genomic bed file")
    parser.add_argument("--positive_bed", required=True, help="Path to the positive sequences bed file")
    parser.add_argument("--model_names", nargs='+', required=True, help="List of pretrained model names")
    parser.add_argument("--plot_average", action='store_true', help="Plot average AUPRC of models")
    parser.add_argument("--plot_variance", action='store_true', help="Plot average AUPRC of models with variance")
    parser.add_argument("--output_name", help="Output name for the plot files")
    parser.add_argument("--output_dir", help="Output directory for the plot files")
    parser.add_argument("--calculate_f1", action='store_true', help="Calculate F1 score instead of AUPRC")
    parser.add_argument("--positive_fasta", help="Path to positive sequences FASTA file for threshold calculation")
    parser.add_argument("--kmer_length", type=int, choices=[3, 4, 5, 6], default=3, help="K-mer length for tokenization (3, 4, 5, or 6)")
    args = parser.parse_args()

    # Set default values if output_name or output_dir are not provided
    if not args.output_name:
        args.output_name = os.path.basename(args.genomic_bed).split('_')[0]
    if not args.output_dir:
        args.output_dir = "."

    if args.calculate_f1:
        if not args.positive_fasta:
            raise ValueError("You must provide a FASTA file with positive sequences for threshold calculation.")

    models_data = {
        "1:1": [],
        "1:3": [],
        "shuffled": []
    }

    for model_name in args.model_names:
        tokenizer = AutoTokenizer.from_pretrained(model_name, trust_remote_code=True)
        model = AutoModelForSequenceClassification.from_pretrained(model_name, trust_remote_code=True, num_labels=2).to(device)
        predicted_windows = predict_genomic_regions(args.genomic_bed, args.positive_bed, model, tokenizer, args.kmer_length)

        if args.calculate_f1:
            threshold = calculate_threshold(args.positive_fasta, model, tokenizer, args.kmer_length)
            f1 = calculate_f1_score(args.positive_bed, predicted_windows, threshold)
            print(f"F1 Score for {model_name}: {f1}")
        else:
            precision, recall, auprc = calculate_precision_recall_curve(args.positive_bed, predicted_windows)
            if "Neg1x" in model_name:
                models_data["1:1"].append((model_name, precision, recall, auprc))
            elif "Neg3x" in model_name:
                models_data["1:3"].append((model_name, precision, recall, auprc))
            else:
                models_data["shuffled"].append((model_name, precision, recall, auprc))

    if not args.calculate_f1:
        if args.plot_average:
            plot_average_auprc(models_data, args.output_name, args.output_dir)
        if args.plot_variance:
            plot_variance_auprc(models_data, args.output_name, args.output_dir)
        else:
            # Plot individual curves
            plot_precision_recall_curves(models_data["1:1"], f"{args.output_name}_1to1", args.output_dir)
            plot_precision_recall_curves(models_data["1:3"], f"{args.output_name}_1to3", args.output_dir)
            plot_precision_recall_curves(models_data["shuffled"], f"{args.output_name}_shuffled", args.output_dir)

if __name__ == "__main__":
    main()

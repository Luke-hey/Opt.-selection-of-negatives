import argparse
from transformers import AutoTokenizer, AutoModelForSequenceClassification
import torch
from sklearn.metrics import f1_score
from Bio import SeqIO

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")


def load_fasta_sequences(file_path):
    sequences = []
    with open(file_path, 'r') as fasta_file:
        for record in SeqIO.parse(fasta_file, 'fasta'):
            sequence = str(record.seq).replace('U', 'T')
            sequence = sequence.upper()
            # Insert a space after every 3 characters
            spaced_sequence = ' '.join([sequence[i:i + 3] for i in range(0, len(sequence), 3)])
            sequences.append(spaced_sequence)
    return sequences


def load_sequences_with_labels(positive_file, negative_file):
    positive_sequences = load_fasta_sequences(positive_file)
    negative_sequences = load_fasta_sequences(negative_file)
    return positive_sequences, negative_sequences


def pred(seq, model, tokenizer, threshold=0.90):
    with torch.no_grad():
        model.eval()
        inputs = tokenizer(seq, return_tensors="pt").to(device)
        outputs = model(**inputs)
        logits = outputs.logits
        probabilities = torch.nn.functional.softmax(logits, dim=-1)
        prediction = (probabilities[:, 1] > threshold).long()
        return prediction.item()


def main():
    parser = argparse.ArgumentParser(description="DNA Sequence Classification Script",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--model_name", required=True, help="Fine-tuned model name in Hugging Face Model Hub")
    parser.add_argument("--positive_file", required=True, help="Input positive sequence file in FASTA format")
    parser.add_argument("--negative_file", required=True, help="Input negative sequence file in FASTA format")
    args = parser.parse_args()

    model_name = args.model_name
    positive_file = args.positive_file
    negative_file = args.negative_file

    tokenizer = AutoTokenizer.from_pretrained(model_name, trust_remote_code=True)
    model = AutoModelForSequenceClassification.from_pretrained(model_name, trust_remote_code=True, num_labels=2).to(device)

    positive_sequences, negative_sequences = load_sequences_with_labels(positive_file, negative_file)

    true_labels = [1] * len(positive_sequences) + [0] * len(negative_sequences)
    predicted_labels = []

    all_sequences = positive_sequences + negative_sequences
    for sequence in all_sequences:
        prediction = pred(sequence, model, tokenizer)
        predicted_labels.append(prediction)

    f1 = f1_score(true_labels, predicted_labels)

    print(f"F1 Score: {f1}")


if __name__ == "__main__":
    main()

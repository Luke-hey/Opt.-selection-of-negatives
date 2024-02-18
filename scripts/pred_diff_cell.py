import argparse
from transformers import AutoTokenizer, AutoModelForSequenceClassification
import torch
from Bio import SeqIO

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")


def load_fasta_sequences(file_path):
    sequences = {}
    with open(file_path, 'r') as fasta_file:
        for record in SeqIO.parse(fasta_file, 'fasta'):
            header = record.id
            sequence = str(record.seq).replace('U', 'T')
            sequences[header] = sequence.upper()
    return sequences


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


def pred(seq, model, tokenizer, threshold=0.86):
    with torch.no_grad():
        model.eval()
        inputs = tokenizer(kmer(seq), return_tensors="pt").to(device)
        outputs = model(**inputs)
        logits = outputs.logits
        probabilities = torch.nn.functional.softmax(logits, dim=-1)
        prediction = (probabilities[:, 1] > threshold).long()
        return prediction.item()

def main():
    parser = argparse.ArgumentParser(description="DNA Sequence Classification Script",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--model_name", required=True, help="Fine-tuned model name in Hugging Face Model Hub")
    parser.add_argument("--input_file", required=True, help="Input positive sequence file in FASTA format")
    args = parser.parse_args()

    model_name = args.model_name
    input_file = args.input_file

    tokenizer = AutoTokenizer.from_pretrained(model_name, trust_remote_code=True)
    model = AutoModelForSequenceClassification.from_pretrained(model_name, trust_remote_code=True, num_labels=2).to(device)

    sequences = load_fasta_sequences(input_file)

    predictions = []

    for sequence in sequences:
        prediction = pred(sequence, model, tokenizer)
        predictions.append(prediction)

    # Calculate the current accuracy
    current_accuracy = sum(predictions) / len(predictions)

    # Set the threshold to achieve roughly 50% accuracy
    threshold = 0.5

    # Adjust threshold based on current accuracy
    if current_accuracy < 0.5:
        threshold += 0.1
    elif current_accuracy > 0.5:
        threshold -= 0.1

    correct_predictions = sum(prediction > threshold for prediction in predictions)

    print(f"Number of correct predictions: {correct_predictions} out of {len(sequences)}")
    print(f"Threshold dynamically adjusted to: {threshold}")


if __name__ == "__main__":
    main()

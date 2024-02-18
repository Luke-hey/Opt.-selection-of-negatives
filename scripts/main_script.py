import subprocess
import argparse
import os

def run_pad_positives(input_pos, flank_size):
    # Run pad_positives.py with the specified input and flank size
    command = ['python', 'pad_positives.py', input_pos, '--flank_size', str(flank_size)]
    subprocess.run(command)

def run_bed_to_fasta(input_pos):
    # Run bed_to_fasta.py with the specified input
    output_fasta = os.path.splitext(os.path.basename(input_pos))[0] + "_adjusted.fa"
    command = ['python', 'bed_to_fasta.py', input_pos, '--output_file', output_fasta]
    subprocess.run(command)

def run_bedtools_script(input_pos, input_low_conf1, input_low_conf2, output_prefix):
    # Run bedtools_script.py with the specified inputs and output prefix
    command = ['python', 'bedtools_script.py', input_pos, input_low_conf1, input_low_conf2, output_prefix]
    subprocess.run(command)

def run_sampling_script(genes_without_pos_file, genes_without_pos, positives_file, num_samples, output_prefix, sampling_mode):
    # Run the appropriate sampling script based on the sampling mode
    if sampling_mode == 'hard':
        output_file = f"{output_prefix}_hard_negatives_output.fa"
        command = ['python', 'sampling_hardneg.py', genes_without_pos, positives_file, output_file, '--num_sequences_to_save', str(num_samples)]
        subprocess.run(command)
    else:
        output_file = f"{output_prefix}_normal_negatives_output.bed"
        command = ['python', 'generate_negatives.py', genes_without_pos_file, output_file, '--num_samples', str(num_samples)]
        subprocess.run(command)

def main():
    parser = argparse.ArgumentParser(description="Run a series of scripts for processing biological sequence data.")
    parser.add_argument("--input_pos", required=True, help="Path to the input BED file for positive sequences.")
    parser.add_argument("--input_low_conf1", required=True, help="Path to the input BED file for low confidence sequences 1.")
    parser.add_argument("--input_low_conf2", required=True, help="Path to the input BED file for low confidence sequences 2.")
    parser.add_argument("--num_samples", type=int, required=True, help="Number of sequences to be sampled.")
    parser.add_argument("--output_prefix", required=True, help="Prefix for resulting files.")
    parser.add_argument("--sampling_mode", choices=['hard', 'normal'], required=True, help="Sampling mode: 'hard' for hard negatives or 'normal' for normal negatives.")

    args = parser.parse_args()

    # Step 1: Run pad_positives.py
    run_pad_positives(args.input_pos, flank_size=50)

    # Step 2: Run bed_to_fasta.py
    run_bed_to_fasta(args.input_pos)

    # Step 3: Run bedtools_script.py
    run_bedtools_script(args.input_pos, args.input_low_conf1, args.input_low_conf2, args.output_prefix)
    run_bed_to_fasta(f"{args.output_prefix}_genes_without_pos.bed")
    # Step 4: Run sampling script (hard or normal)
    run_sampling_script(genes_without_pos_file=f"{args.output_prefix}_genes_without_pos.bed", genes_without_pos=f"{args.output_prefix}_genes_without_pos.fa", positives_file=f"{args.output_prefix}_high_conf_pos.fa", num_samples=args.num_samples, output_prefix=args.output_prefix, sampling_mode=args.sampling_mode)

if __name__ == "__main__":
    main()

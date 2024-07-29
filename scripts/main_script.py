import subprocess
import argparse
import os

def get_script_path(script_name):
    # Get the directory of the currently executing script
    current_dir = os.path.dirname(os.path.realpath(__file__))
    # Construct the path to the script
    return os.path.join(current_dir, script_name)

def count_positive_sequences(input_pos):
    # Count the number of positive sequences in the BED file
    with open(input_pos, 'r') as file:
        count = sum(1 for line in file)
    return count

def run_pad_positives(input_pos):
    # Run pad_positives.py with the specified input and flank size
    command = ['python', get_script_path('new_padd_pos.py'), input_pos]
    subprocess.run(command)

def run_bed_to_fasta(input_pos):
    # Run bed_to_fasta.py with the specified input
    output_fasta = os.path.splitext(os.path.basename(input_pos))[0] + "_adjusted.fa"
    command = ['python', get_script_path('bed_to_fasta.py'), input_pos, '--output_file', output_fasta]
    subprocess.run(command)

def run_bedtools_script(input_pos, input_low_conf1, input_low_conf2, test_sites, test_regions, output_prefix):
    # Run bedtools_script.py with the specified inputs and output prefix
    command = ['python', get_script_path('bedtools_script.py'), input_pos, input_low_conf1, input_low_conf2, test_sites, test_regions, output_prefix]
    subprocess.run(command)

def setup_processing(input_pos, input_low_conf1, input_low_conf2, test_sites, test_regions, output_prefix):
    # Check if the preprocessing steps have already been done
    # If not, perform the preprocessing steps
    genes_without_pos_bed = f"{output_prefix}_genes_without_pos.bed"
    #genes_without_pos_fa = f"{output_prefix}_genes_without_pos.fa"

    if not (os.path.exists(genes_without_pos_bed)):
        run_pad_positives(input_pos)
        #run_bed_to_fasta(genes_without_pos_bed)
    if not (os.path.exists(f"{output_prefix}_train_pos.bed")):
        run_bedtools_script(input_pos, input_low_conf1, input_low_conf2, test_sites, test_regions, output_prefix)
        run_bed_to_fasta(f"{output_prefix}_train_pos.bed")

def run_sampling_script(genes_without_pos_file, genes_without_pos, positives_file, positives_bed, input_low_conf1, input_low_conf2, num_samples, existing_negatives, output_prefix, sampling_mode, mult):
    # Run the appropriate sampling script based on the sampling mode
    if sampling_mode == 'hard':
        output_file = f"{output_prefix}_hard_negatives_output{num_samples}.fa"
        command = ['python', get_script_path('sampling_hardneg.py'), genes_without_pos, positives_file, output_file, '--num_sequences_to_save', str(num_samples)]
        subprocess.run(command)
    else:
        output_file = f"{output_prefix}_Neg{mult}x.bed"
        command = ['python', get_script_path('generate_negatives.py'), genes_without_pos_file, positives_bed, input_low_conf1, input_low_conf2, output_file, '--num_samples', str(num_samples)]
        if existing_negatives:
            command += ['--existing_negatives'] + existing_negatives
        subprocess.run(command)
        run_bed_to_fasta(output_file)

def main():
    parser = argparse.ArgumentParser(description="Run a series of scripts for processing biological sequence data.")
    parser.add_argument("--input_pos", required=True, help="Path to the input BED file for positive sequences.")
    parser.add_argument("--input_low_conf1", required=True, help="Path to the input BED file for low confidence sequences 1.")
    parser.add_argument("--input_low_conf2", required=True, help="Path to the input BED file for low confidence sequences 2.")
    parser.add_argument("--test_sites", required=True, help="Path to the input BED file for test sites.")
    parser.add_argument("--test_regions", required=True, help="Path to the input BED file for test regions.")
    parser.add_argument("--output_prefix", required=True, help="Prefix for resulting files.")
    parser.add_argument("--sampling_mode", choices=['hard', 'normal'], required=True, help="Sampling mode: 'hard' for hard negatives or 'normal' for normal negatives.")
    parser.add_argument("--num_samples", type=int, help="Multiplyer for number of negatives to be sampled. If not provided, defaults to the number of positive sequences in the input BED file.")
    parser.add_argument("--existing_negatives", nargs='+', default=[], help="List of existing negative files to avoid overlaps")

    args = parser.parse_args()

    setup_processing(args.input_pos, args.input_low_conf1, args.input_low_conf2, args.test_sites, args.test_regions, args.output_prefix)

    num_samples = args.num_samples * count_positive_sequences(f"{args.output_prefix}_train_pos.bed")
    if num_samples is None:
        num_samples = count_positive_sequences(f"{args.output_prefix}_train_pos.bed")

    # Step 4: Run sampling script (hard or normal)
    run_sampling_script(genes_without_pos_file=f"{args.output_prefix}_genes_without_pos.bed", genes_without_pos=f"{args.output_prefix}_genes_without_pos.fa", positives_file = f"{args.input_pos.split('.')[0]}.fa", positives_bed = f"{args.input_pos.split('.')[0]}.bed", input_low_conf1=args.input_low_conf1, input_low_conf2=args.input_low_conf2, num_samples=num_samples, existing_negatives=args.existing_negatives, output_prefix=args.output_prefix, sampling_mode=args.sampling_mode, mult=args.num_samples)

if __name__ == "__main__":
    main()

#!/usr/bin/env python

import argparse
import pybedtools
from pybedtools import BedTool
import os


# Get the directory of the current script
script_dir = os.path.dirname(os.path.abspath(__file__))
# Construct the relative path
REFERENCE_GENOME = os.path.join(script_dir, "../hg38.fa")
#REFERENCE_GENOME = "/home/bubu23/hg38.fa"        # add path here

def extract_sequences_from_bed(bed_file, output_file):
    bed_file = BedTool(bed_file)
    hg_fasta = BedTool(REFERENCE_GENOME)
    sequences = bed_file.sequence(fi=hg_fasta, s=True, nameOnly=True)

    # Open the generated FASTA file for reading
    with open(sequences.seqfn) as fasta_file:
        # Open the output file for writing
        with open(output_file, 'w') as output_handle:
            sequence = ''
            for line in fasta_file:
                # If it's a sequence header line
                if line.startswith('>'):
                    # Write the previous sequence to the output file in uppercase
                    if sequence:
                        output_handle.write(sequence.upper() + '\n')
                        sequence = ''
                    output_handle.write(line)
                else:
                    sequence += line.strip()
            # Write the last sequence to the output file in uppercase
            if sequence:
                output_handle.write(sequence.upper() + '\n')

    # Print the number of sequences in the resulting FASTA file
    num_sequences = sum(1 for _ in open(output_file)) // 2
    print(f"Number of sequences in {output_file}: {num_sequences}")

    return output_file


def main():
    parser = argparse.ArgumentParser(description="Extract sequences from BED file using a fixed reference genome.")
    parser.add_argument("bed_file", help="Input BED file")
    parser.add_argument("--output_file", "-o", default="output.fa", help="Output file name (default: output.fa)")

    args = parser.parse_args()

    # Remove the old extension and append .fa to the input BED file name for the output file
    input_basename = os.path.splitext(os.path.basename(args.bed_file))[0]
    output_file = input_basename + ".fa"

    extract_sequences_from_bed(args.bed_file, output_file)

if __name__ == "__main__":
    main()

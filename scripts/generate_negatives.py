#!/usr/bin/env python

import argparse
import pybedtools
from pybedtools import BedTool
import random

SEQUENCE_LENGTH = 101  # Set your desired sequence length

def sample_non_overlapping_negative_sequences(input_bed, output_bed, num_samples=12558):
    # Convert the input BED file to a BedTool object
    input_bedtool = BedTool(input_bed)

    # Get a list of intervals from the input BED file
    input_intervals = list(input_bedtool)

    sampled_negative_intervals_bedtool = BedTool([])
    sampled_intervals = set()
    count = 1

    for i in range(num_samples):
        # Shuffle the list of intervals
        random.shuffle(input_intervals)

        # Iterate through the shuffled intervals
        for interval in input_intervals:
            # Get the start, end, chromosome, and strand information of the interval
            start = interval.start
            end = interval.end
            chrom = interval.chrom
            strand = interval.strand if interval.strand else '.'  # Handle cases where strand information is missing

            # Calculate the midpoint of the interval
            midpoint = (start + end) // 2

            # Adjust the interval length to the constant sequence length
            adjusted_start = midpoint - SEQUENCE_LENGTH // 2
            adjusted_end = midpoint + (SEQUENCE_LENGTH // 2) + (SEQUENCE_LENGTH % 2)

            # Check for overlap with already sampled intervals
            if not any((chrom == e[0] and adjusted_start < e[2] and adjusted_end > e[1] and strand == e[3]) for e in sampled_intervals):
                # Preserve the additional columns (4th to last) in the adjusted interval
                adjusted_interval = f"{chrom}\t{adjusted_start}\t{adjusted_end}\t rand_neg_{count}\t{interval[4]}\t{interval[5]}\n"
                count += 1

                # Add the adjusted interval to the output BedTool and update the sampled intervals set
                sampled_negative_intervals_bedtool = sampled_negative_intervals_bedtool.cat(BedTool(adjusted_interval, from_string=True), postmerge=False)
                sampled_intervals.add((chrom, adjusted_start, adjusted_end, strand))
                break

    # Save the sampled negative intervals to the output BED file
    sampled_negative_intervals_bedtool.saveas(output_bed)

def main():
    parser = argparse.ArgumentParser(description="Sample non-overlapping negative sequences from an input BED file.")
    parser.add_argument("input_bed", help="Input BED file")
    parser.add_argument("output_bed", help="Output BED file")
    parser.add_argument("--num_samples", type=int, default=12558, help="Number of sequences to sample (default: 12558)")

    args = parser.parse_args()

    sample_non_overlapping_negative_sequences(args.input_bed, args.output_bed, args.num_samples)

if __name__ == "__main__":
    main()

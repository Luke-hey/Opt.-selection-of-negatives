#!/usr/bin/env python

import argparse
import pybedtools
from pybedtools import BedTool
import random

SEQUENCE_LENGTH = 101  # Set your desired sequence length

def sample_non_overlapping_negative_sequences(input_bed, positive_bed, low_conf_pos, output_bed, num_samples=None, existing_negatives=[]):
    if num_samples is None:
        num_samples = sum(1 for _ in open(positive_bed))
    # Convert the input BED files to BedTool objects
    input_bedtool = BedTool(input_bed)
    positive_bedtool = BedTool(positive_bed)

    # Merge existing negative files into a single BedTool object
    existing_negatives_bed_string = ''
    for neg_file in existing_negatives:
        with open(neg_file, 'r') as f:
            existing_negatives_bed_string += f.read()

    # Create BedTool object from the concatenated string
    existing_negatives_bedtool = BedTool(existing_negatives_bed_string, from_string=True)

    sampled_negative_intervals_bedtool = BedTool([])
    sampled_intervals = set()

    while len(sampled_intervals) < num_samples:
        # Shuffle the list of intervals from the input BED file
        input_intervals = list(input_bedtool)
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

            # Check for overlap with positive intervals and existing negative intervals
            if not (positive_bedtool.intersect(BedTool(f"{chrom}\t{adjusted_start}\t{adjusted_end}\t{interval[3]}\t{interval[4]}\t{interval[5]}\n", from_string=True), u=True, s=True) or
                    low_conf_pos.intersect(BedTool(f"{chrom}\t{adjusted_start}\t{adjusted_end}\t{interval[3]}\t{interval[4]}\t{interval[5]}\n", from_string=True), u=True, s=True) or
                    existing_negatives_bedtool.intersect(BedTool(f"{chrom}\t{adjusted_start}\t{adjusted_end}\t{interval[3]}\t{interval[4]}\t{interval[5]}\n", from_string=True), u=True, s=True)):
                # Check for overlap with already sampled negative intervals
                overlap = False
                for sampled_interval in sampled_intervals:
                    if chrom == sampled_interval[0] and \
                       max(adjusted_start, sampled_interval[1]) < min(adjusted_end, sampled_interval[2]):
                        overlap = True
                        break

                # If no overlap, add the adjusted interval to the output BedTool and update the sampled intervals set
                if not overlap:
                    sampled_intervals.add((chrom, adjusted_start, adjusted_end, strand))
                    sampled_negative_intervals_bedtool = sampled_negative_intervals_bedtool.cat(BedTool(f"{chrom}\t{adjusted_start}\t{adjusted_end}\t{interval[3]}\t{interval[4]}\t{interval[5]}\n", from_string=True), postmerge=False)

            if len(sampled_intervals) >= num_samples:
                break

    # Save the sampled negative intervals to the output BED file
    sampled_negative_intervals_bedtool.saveas(output_bed)

def main():
    parser = argparse.ArgumentParser(description="Sample non-overlapping negative sequences from an input BED file.")
    parser.add_argument("input_bed", help="Input BED file containing regions where positive sequences have been subtracted out")
    parser.add_argument("positive_bed", help="BED file containing positive sequences")
    parser.add_argument("pos_low_conf1", help="BED file containing low positive sequences")
    parser.add_argument("pos_low_conf2", help="BED file containing low positive sequences")
    parser.add_argument("output_bed", help="Output BED file")
    parser.add_argument("--num_samples", type=int, default=None, help="Number of sequences to sample (default: number of positive sequences in positive_bed)")
    parser.add_argument("--existing_negatives", nargs='+', default=[], help="List of existing negative files to avoid overlaps")

    args = parser.parse_args()
    low_conf_pos = BedTool(args.pos_low_conf1).cat(BedTool(args.pos_low_conf2), postmerge=False)
    sample_non_overlapping_negative_sequences(args.input_bed, args.positive_bed, low_conf_pos, args.output_bed, args.num_samples, args.existing_negatives)

if __name__ == "__main__":
    main()

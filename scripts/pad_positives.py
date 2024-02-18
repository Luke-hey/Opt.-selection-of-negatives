#!/usr/bin/env python
from pybedtools import BedTool
import argparse
import os

count = 0
def adjust_interval(interval, flank_size):
    global count
    try:
        # Calculate the center of the interval
        center = (interval.start + interval.end) // 2
        count += 1
        # Adjust start and end positions based on strand
        adjusted_start = max(0, center - flank_size)
        adjusted_end = center + flank_size + 1


        adjusted_name = f"positive_{count + 1}"

        # Create a new formatted string for the adjusted interval
        adjusted_interval_str = f"{interval.chrom}\t{int(adjusted_start)}\t{int(adjusted_end)}\tpositive_{count}\t{interval.score}\t{interval.strand}"

        # Create a string for the original interval
        original_interval_str = str(interval)

        return f"{adjusted_interval_str}" #"\t{original_interval_str}"

    except Exception as e:
        print(f"Error adjusting interval: {interval}, Error: {e}")
        raise

def expand_bed_file(input_bed, flank_size=50):
    # Load and sort the input BED file
    bed = BedTool(input_bed).sort()

    # Adjust intervals without parallel processing
    adjusted_intervals = [adjust_interval(interval, flank_size) for interval in bed]

    # Create a new BED file with adjusted intervals and original ranges
    adjusted_bed = BedTool("\n".join(adjusted_intervals), from_string=True)

    # Generate the output BED file name
    output_bed_file = os.path.basename(input_bed)

    # Save the adjusted BED file
    adjusted_bed.saveas(output_bed_file)
    print(f"Adjusted BED file saved as {output_bed_file}")

def parse_args():
    parser = argparse.ArgumentParser(description="Adjust intervals in a BED file by adding flanks.")
    parser.add_argument("input_bed", help="Path to the input BED file.")
    parser.add_argument("--flank_size", type=int, default=50, help="Size of flanks to add on both sides of the intervals.")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()
    input_bed_file = args.input_bed
    flank_size = args.flank_size

    expand_bed_file(input_bed_file, flank_size)

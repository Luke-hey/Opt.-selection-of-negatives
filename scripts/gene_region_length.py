import argparse
import pybedtools

def calculate_total_gene_region_length(bed_file):
    # Load the BED file using pybedtools
    bed = pybedtools.BedTool(bed_file)

    # Initialize total length counter
    total_length = 0

    # Calculate the length of each gene region and accumulate the total length
    for interval in bed:
        gene_length = len(interval)
        total_length += gene_length

    # Print out the total length of gene regions
    print(f"Total length of gene regions: {total_length} bp")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate total length of gene regions from a BED file")
    parser.add_argument("bed_file", help="Path to the BED file")
    args = parser.parse_args()

    # Call the function with the provided BED file
    calculate_total_gene_region_length(args.bed_file)

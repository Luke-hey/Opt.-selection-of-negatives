#!/usr/bin/env python

import argparse
import pybedtools
from pybedtools import BedTool

ANNOTATIONS_FILE = "/home/bubu23/gencode.v44.annotation.gtf.gz"  # add path to annotation here
HG_FASTA_FILE = "/home/bubu23/hg38.fa"                           # add path to genome here

def process_bed_files(positives_file, low_conf_pos1_file, low_conf_pos2_file, test_sites_file, test_regions_file, output_prefix):
    # Read constant bed files
    annotations = BedTool(ANNOTATIONS_FILE)
    hg_fasta = BedTool(HG_FASTA_FILE)

    # Read variable bed files
    positives = BedTool(positives_file)
    low_conf_pos1 = BedTool(low_conf_pos1_file)
    low_conf_pos2 = BedTool(low_conf_pos2_file)
    test_sites = BedTool(test_sites_file)
    test_regions = BedTool(test_regions_file)


    train_positives = positives.intersect(test_sites, s=True, v=True)
    train_positives.saveas(f"{output_prefix}_train_pos.bed")
    # Filter rows based on the specified feature type
    feature_type_to_filter = "gene"
    filtered_genes = annotations.filter(lambda x: x[2] == feature_type_to_filter)

    # Replace the feature type in the third column with the value in the last column
    modified_genes = [
        [x[0], x[3], x[4], x[8].split(";")[0][8:], x[5], x[6]]
        for x in filtered_genes
    ]

    # Create a new BedTool with the modified information (gene regions on genome)
    modified_genes_bed = BedTool(modified_genes)

    # Intersect positives with the gene regions
    genes_with_pos = modified_genes_bed.intersect(positives, u=True, s=True)
    #genes_with_pos.saveas(f"{output_prefix}_genes_with_pos.bed")

    # Concatenate low_conf_pos1 and low_conf_pos2
    low_conf_pos = low_conf_pos1.cat(low_conf_pos2, postmerge=False)
    #low_conf_pos.saveas(f"{output_prefix}_low_conf_pos.bed")

    # Subtract positives and low_conf_pos from positive gene regions
    #genes_without_pos = modified_genes_bed.subtract(positives, s=True)
    genes_without_pos = genes_with_pos.subtract(positives, s=True)
    genes_without_pos = genes_without_pos.subtract(low_conf_pos, s=True)
    genes_without_pos = genes_without_pos.subtract(test_regions, s=True)
    genes_without_pos = genes_without_pos.subtract(test_sites, s=True)
    genes_without_pos.saveas(f"{output_prefix}_genes_without_pos.bed")


def main():
    parser = argparse.ArgumentParser(description="Process BED files and generate output files.")
    parser.add_argument("positives_file", help="Input BED file for positives")
    parser.add_argument("low_conf_pos1_file", help="Input BED file for low_conf_pos1")
    parser.add_argument("low_conf_pos2_file", help="Input BED file for low_conf_pos2")
    parser.add_argument("testing_sites", help="Input BED file for testing sites")
    parser.add_argument("testing_regions", help="Input BED file for testing regions")
    parser.add_argument("prefix", help="prefix of output files")

    args = parser.parse_args()

    # Extract output prefix from the positive file until the second underscore
     # Extract output prefix from the positive file using the file name without extension
    #output_prefix = os.path.splitext(os.path.basename(args.positives_file))[0]

    process_bed_files(args.positives_file, args.low_conf_pos1_file, args.low_conf_pos2_file, args.testing_sites, args.testing_regions, args.prefix)

if __name__ == "__main__":
    main()

# generate fasta file
#modified_genes_bed = modified_genes_bed.sequence(fi=hg_fasta, s = True, name = True, fo="genes.fa")
#print(open(modified_genes_bed.seqfn).read())

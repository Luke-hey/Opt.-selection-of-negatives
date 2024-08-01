#!/bin/bash

# Check if at least one protein name is provided as an argument
if [ "$#" -lt 1 ]; then
    echo "Usage: $0 <protein_name> [<protein_name> ...]"
    exit 1
fi

# Variables
script="../scripts/DNABERT_3mer.py"
directories=("Run_1" "Run_2" "Run_3")
queue_file="protein_queue.txt"

# Function to process directories
process_directories() {
    local protein="$1"
    local protein_dir="./$protein"
    local pos_fasta="${protein_dir}/${protein}_train_pos.fa"
    for dir in "${directories[@]}"; do
        # Find the negative FASTA files in the current directory
        neg1x_fasta="${protein_dir}/${dir}/${protein}_Neg1x.fa"
        neg3x_fasta="${protein_dir}/${dir}/${protein}_Neg3x.fa"
        shuffled_fasta="${protein_dir}/${dir}/${protein}_shuffled.fa"

        # Define output directories
        output_dir="${protein_dir}/${dir}/"

        # Check and process Neg1x FASTA file
        if [[ -f "$neg1x_fasta" ]]; then
            if [ ! -d "${output_dir}/finetuned_DNABERT3mer_${protein}_Neg1x" ]; then
                echo "Processing $neg1x_fasta"
                if ! python "$script" --positive_sequences "$pos_fasta" --negative_sequences "$neg1x_fasta" --output_dir "$output_dir"; then
                    echo "Error processing $neg1x_fasta"
                fi
            else
                echo "Output directory ${output_dir}/finetuned_DNABERT3mer_${protein}_Neg1x already exists. Skipping."
            fi
        else
            echo "Negative FASTA file $neg1x_fasta not found"
        fi

        # Check and process Neg3x FASTA file
        if [[ -f "$neg3x_fasta" ]]; then
            if [ ! -d "${output_dir}/finetuned_DNABERT3mer_${protein}_Neg3x" ]; then
                echo "Processing $neg3x_fasta"
                if ! python "$script" --positive_sequences "$pos_fasta" --negative_sequences "$neg3x_fasta" --output_dir "$output_dir"; then
                    echo "Error processing $neg3x_fasta"
                fi
            else
                echo "Output directory ${output_dir}/finetuned_DNABERT3mer_${protein}_Neg3x already exists. Skipping."
            fi
        else
            echo "Negative FASTA file $neg3x_fasta not found"
        fi

        # Check and process Shuffled FASTA file
        if [[ -f "$shuffled_fasta" ]]; then
            if [ ! -d "${output_dir}/finetuned_DNABERT3mer_${protein}_shuffled" ]; then
                echo "Processing $shuffled_fasta"
                if ! python "$script" --positive_sequences "$pos_fasta" --negative_sequences "$shuffled_fasta" --output_dir "$output_dir"; then
                    echo "Error processing $shuffled_fasta"
                fi
            else
                echo "Output directory ${output_dir}/finetuned_DNABERT3mer_${protein}_shuffled already exists. Skipping."
            fi
        else
            echo "Negative FASTA file $shuffled_fasta not found"
        fi
    done
}

# Create the queue file if it doesn't exist
if [ ! -f "$queue_file" ]; then
    touch "$queue_file"
fi

# Add each provided protein name to the queue if it's not already there
for protein_name in "$@"; do
    if ! grep -q "^${protein_name}$" "$queue_file"; then
        echo "$protein_name" >> "$queue_file"
    fi
done

# Process the queue
while IFS= read -r protein; do
    echo "Processing protein: $protein"
    process_directories "$protein"
    # Remove the processed protein from the queue
    sed -i.bak "/^${protein}$/d" "$queue_file" && rm "$queue_file.bak"
done < "$queue_file"

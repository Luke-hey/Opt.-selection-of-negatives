#!/bin/bash

# Usage function to display help
usage() {
    echo "Usage: $0 -p <protein_name> -c <cell_type> -n <num_samples>"
    echo "  -p <protein_name>   Name of the protein"
    echo "  -c <cell_type>      Type of the cell"
    echo "  -n <num_samples>    Multiplier for generating negatives. For example, if n=2, the number of negatives will be twice the number of positives."
    exit 1
}

# Parse command line arguments
while getopts ":p:c:n:" opt; do
  case $opt in
    p) PROTEIN_NAME="$OPTARG"
    ;;
    c) CELL_TYPE="$OPTARG"
    ;;
    n) NUM_SAMPLES="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
        usage
    ;;
  esac
done

# Check if all arguments are provided
if [ -z "$PROTEIN_NAME" ] || [ -z "$CELL_TYPE" ] || [ -z "$NUM_SAMPLES" ]; then
    usage
fi

# File paths and common prefix
INPUT_POS="${PROTEIN_NAME}_${CELL_TYPE}_1_2.bed"
INPUT_LOW_CONF1="${PROTEIN_NAME}_${CELL_TYPE}_1.bed"
INPUT_LOW_CONF2="${PROTEIN_NAME}_${CELL_TYPE}_2.bed"
TEST_SITES="${PROTEIN_NAME}_sites.bed"
TEST_REGIONS="${PROTEIN_NAME}_regions.bed"
OUTPUT_PREFIX="${PROTEIN_NAME}"

# Ensure directories for runs exist
for RUN_DIR in Run_1 Run_2 Run_3; do
    if [ ! -d "$RUN_DIR" ]; then
        mkdir -p "$RUN_DIR"
    fi
done

# Function to check if negative files already exist
negative_files_exist() {
    local run_dir=$1
    local num_samples=$2
    local neg_fasta="${run_dir}/${PROTEIN_NAME}_Neg${num_samples}x.fa"
    local neg_bed="${run_dir}/${PROTEIN_NAME}_Neg${num_samples}x.bed"

    if [[ -f "$neg_fasta" && -f "$neg_bed" ]]; then
        return 0
    else
        return 1
    fi
}

# Run the commands in a loop
for RUN in {1..3}; do
    RUN_DIR="Run_$RUN"
    echo "Running iteration $RUN in directory $RUN_DIR"

    # Check if the negative files already exist in the current run directory
    if negative_files_exist "$RUN_DIR" "$NUM_SAMPLES"; then
        echo "Negative files already exist in $RUN_DIR. Skipping iteration $RUN."
        continue
    fi

    # Print current directory and list files
    echo "Current directory: $(pwd)"
    echo "Files in current directory before running the script:"
    ls -l

    # Use relative path for the script
    python ./scripts/main_script.py --input_pos "$INPUT_POS" --input_low_conf1 "$INPUT_LOW_CONF1" --input_low_conf2 "$INPUT_LOW_CONF2" --test_sites "$TEST_SITES" --test_regions "$TEST_REGIONS" --output_prefix "$OUTPUT_PREFIX" --sampling_mode normal --num_samples "$NUM_SAMPLES"

    if [ $? -ne 0 ]; then
        echo "Error: Python script failed during iteration $RUN"
        exit 1
    fi

    # Add a delay to ensure files are fully generated
    sleep 2

    # Print current directory and list files after running the script
    echo "Files in current directory after running the script:"
    ls -l

    # Verify the files to be moved exist
    files_to_move=(${OUTPUT_PREFIX}_Neg*)
    if [ ${#files_to_move[@]} -eq 0 ]; then
        echo "Error: No files found to move for pattern ${OUTPUT_PREFIX}_Neg*"
        exit 1
    fi

    # Move the output files
    for file in "${files_to_move[@]}"; do
        mv "$file" "$RUN_DIR/"
        if [ $? -ne 0 ]; then
            echo "Error: Failed to move file $file during iteration $RUN"
            exit 1
        fi
    done
done

echo "Script execution completed"

#!/bin/bash

# Usage function to display help
usage() {
    echo "Usage: $0 <protein_name> [<protein_name> ...]"
    exit 1
}

# Check if at least one protein name is provided as an argument
if [ "$#" -lt 1 ]; then
    usage
fi

# Variables
script_dir=$(dirname "$0")
script="$script_dir/pred_genomic_region.py"
directories=("Run_1" "Run_2" "Run_3")
model_suffixes=("Neg1x" "Neg3x" "shuffled")
queue_file="plot_queue.txt"

# Function to process directories for plotting
process_directories() {
    local protein="$1"
    local protein_dir="./$protein"
    local genomic_bed="${protein_dir}/${protein}_regions.bed"
    local positive_bed="${protein_dir}/${protein}_sites.bed"
    local model_names=()

    echo "Processing directories for protein: $protein"
    echo "Checking existence of $genomic_bed and $positive_bed"

    # Check if genomic_bed and positive_bed files exist
    if [[ ! -f "$genomic_bed" ]] || [[ ! -f "$positive_bed" ]]; then
        echo "Genomic or positive BED file not found for $protein"
        return 1
    fi

    for dir in "${directories[@]}"; do
        for suffix in "${model_suffixes[@]}"; do
            model_dir="${protein_dir}/${dir}/finetuned_DNABERT3mer_${protein}_${suffix}"
            if [ -d "$model_dir" ]; then
                echo "Found model directory: $model_dir"
                model_names+=("$model_dir")
            else
                echo "Model directory not found: $model_dir"
            fi
        done
    done

    # Check if model_names array is empty
    if [ ${#model_names[@]} -eq 0 ]; then
        echo "No model directories found for $protein"
        return 1
    fi

    # Call the plot_script.py script with the necessary arguments
    python $script --genomic_bed "$genomic_bed" --positive_bed "$positive_bed" --model_names "${model_names[@]}" --plot_variance --output_name "${protein}_variance_plot" --output_dir "${protein_dir}"
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
    echo "Processing plotting for protein: $protein"
    if process_directories "$protein"; then
        # Remove the processed protein from the queue
        sed -i.bak "/^${protein}$/d" "$queue_file" && rm "$queue_file.bak"
    else
        echo "Failed to process $protein for plotting."
    fi
done < "$queue_file"

#!/bin/bash

# Check if the correct number of arguments are provided
if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <number_of_runs> <protein_dir_name1> [<protein_dir_name2> ...]"
    exit 1
fi

# Variables
number_of_runs="$1"
shift # Remove the first argument, so the rest can be treated as protein directory names

# Create the uber directory named "proteins"
uber_directory_name="proteins"
mkdir -p "$uber_directory_name"

# Create directories for each protein inside the uber directory
for protein_dir_name in "$@"; do
    echo "Setting up ${protein_dir_name} with ${number_of_runs} run directories inside ${uber_directory_name}."
    
    # Create the main protein directory inside the uber directory
    main_protein_dir="${uber_directory_name}/${protein_dir_name}"
    mkdir -p "$main_protein_dir"
    
    # Create the specified number of run directories within the main protein directory
    for ((i = 1; i <= number_of_runs; i++)); do
        run_dir="${main_protein_dir}/Run_${i}"
        mkdir -p "$run_dir"
    done

    echo "Completed: Created ${number_of_runs} run directories in ${main_protein_dir}."
done

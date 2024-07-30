#!/bin/bash

# Check if the correct number of arguments are provided
if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <number_of_runs> <protein_dir_name1> [<protein_dir_name2> ...]"
    exit 1
fi

# Variables
number_of_runs="$1"
shift # Remove the first argument, so the rest can be treated as protein directory names

# Create the uber directory named "proteins" if it doesn't exist
uber_directory_name="proteins"
if [ ! -d "$uber_directory_name" ]; then
    mkdir "$uber_directory_name"
    echo "Created directory: $uber_directory_name"
else
    echo "Directory already exists: $uber_directory_name"
fi

# Create directories for each protein inside the uber directory if they don't exist
for protein_dir_name in "$@"; do
    main_protein_dir="${uber_directory_name}/${protein_dir_name}"

    if [ ! -d "$main_protein_dir" ]; then
        echo "Setting up ${protein_dir_name} with ${number_of_runs} run directories inside ${uber_directory_name}."
        mkdir "$main_protein_dir"
        echo "Created directory: $main_protein_dir"
    else
        echo "Directory already exists: $main_protein_dir"
    fi

    # Create the specified number of run directories within the main protein directory if they don't exist
    for ((i = 1; i <= number_of_runs; i++)); do
        run_dir="${main_protein_dir}/Run_${i}"
        if [ ! -d "$run_dir" ]; then
            mkdir "$run_dir"
            echo "Created directory: $run_dir"
        else
            echo "Directory already exists: $run_dir"
        fi
    done

    echo "Completed: Checked/Created ${number_of_runs} run directories in ${main_protein_dir}."
done

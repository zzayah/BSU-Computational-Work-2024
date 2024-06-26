#!/bin/bash

# Get the current directory path
current_directory=$(dirname "$0")

# Loop through all directories in the current directory
for dir in "$current_directory"/*; do
    # Check if the item is a directory
    if [ -d "$dir" ]; then
        # Run the seed_batch_run.sh file in each directory
        if [ -f "$dir/seed_batch_run.sh" ]; then
            cd "$dir"
            sbatch seed_batch_run.sh ./
            cd -
        fi
    fi
done

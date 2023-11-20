#!/bin/bash

# This script was written with the assistance of ChatGPT 3.5
# The script compares the number of reads mapping to a defined region
# to a specified number of random locations of the same size 
# evenly distributed on the chromosome outside of the region.

# Set a fixed random seed (change this value for different runs)
RANDOM_SEED=42
RANDOM=$RANDOM_SEED

# Set the number of repetitions
NUM_REPETITIONS=20

# Path to the tab-delimited text file
# This file contains the path to the indexed .bam-file, the total read counts and position of the dif site
# an example is given in the file bam_statistics.txt
file_path="/home/aap/mnt/Juergen/Projects/vesicles_Ecoli/vesicle_dna/coli/coli_omv_processed/bam_statistics.txt"

# Output file path
output_file="${file_path%/*}/bam_statistics_extended.txt"

# Create a temporary file to store the updated data
temp_file=$(mktemp)

# Function to get chromosome size from BAM file header
get_chromosome_size() {
    local bam_filename="$1"
    local dna_sequence_id="$2"
    local chromosome_size=$(samtools view -H "$bam_filename" | grep -E "^@SQ.*SN:$dna_sequence_id" | awk '{print $3}' | awk -F':' '{print $2}')
    echo "$chromosome_size"
}

# Write the header to the temporary file with lowercase column headers
echo -e "$(head -n 1 "$file_path" | tr '[:upper:]' '[:lower:]')\tcount\tmean\tmedian\tstddev" > "$temp_file"

# Loop through each line in the file (skip the header)
tail -n +2 "$file_path" | while IFS=$'\t' read -r bam_filename dna_sequence_id mapped_reads region_start region_end; do
    # Get the size of the chromosome specified by dna_sequence_id
    chromosome_size=$(get_chromosome_size "$bam_filename" "$dna_sequence_id")

    # Calculate the size of the region
    region_size=$((region_end - region_start))
    
    # Calculate the remaining space outside the specified region
    remaining_space=$((chromosome_size - region_size))
    
    # Initialize an array to store the random counts
    random_counts=()
    
    # Generate random locations and counts
    for ((i=1; i<=NUM_REPETITIONS; i++)); do
        # Divide the remaining space into equal segments
        segment_size=$((remaining_space / NUM_REPETITIONS))
        
        # Calculate the start and end of the segment
        segment_start=$((region_end + 1 + (i - 1) * segment_size))
        segment_end=$((region_end + 1 + i * segment_size))
        
        # Generate a random start position within the segment
        rand_start=$((segment_start + RANDOM % (segment_size + 1)))
        
        # Calculate the corresponding end position
        rand_end=$((rand_start + region_size))
        
        # Use samtools view -c to count reads for the random location
        rand_count=$(samtools view -c "$bam_filename" "$dna_sequence_id:$rand_start-$rand_end")
        
        # Print random location, chromosome ID, and count
        echo "Chromosome ID: $dna_sequence_id, Random Location: $rand_start-$rand_end, Count: $rand_count"
        
        # Append the random count to the array
        random_counts+=("$rand_count")
    done
    
    # Calculate the mean and standard deviation of the random counts
    sum=0
    for rand_count in "${random_counts[@]}"; do
        sum=$((sum + rand_count))
    done
    mean=$((sum / NUM_REPETITIONS))
    
    # Calculate the median
    IFS=$'\n' sorted_counts=($(sort -n <<<"${random_counts[*]}"))
    middle=$((NUM_REPETITIONS / 2))
    if ((NUM_REPETITIONS % 2 == 0)); then
        median=$(( (sorted_counts[middle - 1] + sorted_counts[middle]) / 2 ))
    else
        median="${sorted_counts[middle]}"
    fi
    
    sum_squared_diff=0
    for rand_count in "${random_counts[@]}"; do
        diff=$((rand_count - mean))
        sum_squared_diff=$((sum_squared_diff + diff * diff))
    done
    variance=$((sum_squared_diff / NUM_REPETITIONS))
    stddev=$(bc <<< "scale=2; sqrt($variance)")
    
    # Use samtools view -c to count reads for the specified region
    count=$(samtools view -c "$bam_filename" "$dna_sequence_id:$region_start-$region_end")
    
    # Append the counts, mean, median, and standard deviation to the current line
    echo -e "$bam_filename\t$dna_sequence_id\t$mapped_reads\t$region_start\t$region_end\t$count\t$mean\t$median\t$stddev" >> "$temp_file"
done

# Replace the original file with the temporary file
mv "$temp_file" "$output_file"

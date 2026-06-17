#!/bin/bash

# Global all-atom RMSD calculation script for all proteins.
# Iterates over each protein in config.py, runs cli-global.py with --stat all-atom-rmsd,
# and writes results to CSV and warnings to log files in the results directory.
# Usage: ./cli-global-all-atom-rmsd.sh

data_dir="../data"
results_dir="../results"
logs_dir="../results/logs"

# Create results directory if it doesn't exist
mkdir -p "$results_dir"

# Create logs directory if it doesn't exist
mkdir -p "$logs_dir"

# Array of proteins from config.py
proteins=("PR" "IN" "RT")

echo "Starting global all-atom RMSD calculations for all proteins..."
echo ""

for protein in "${proteins[@]}"; do
    echo "Processing protein: $protein"

    alignment_file="$data_dir/$protein/refs/alignment.fasta"
    protein_lower=$(echo "$protein" | tr '[:upper:]' '[:lower:]')
    output_csv="$results_dir/${protein_lower}-global-all-atom-rmsd.csv"
    output_log="$logs_dir/${protein_lower}-global-all-atom-rmsd.stderr.log"

    # Check if alignment file exists
    if [ ! -f "$alignment_file" ]; then
        echo "  WARNING: Alignment file not found: $alignment_file"
        echo "  Skipping $protein..."
        echo ""
        continue
    fi

    # Run cli-global.py
    echo "  Running cli-global.py..."
    python3 cli-global.py \
        --base-path "$data_dir" \
        --protein "$protein" \
        --stat all-atom-rmsd \
        --alignment "$alignment_file" \
        > "$output_csv" \
        2> "$output_log"

    exit_code=$?

    if [ $exit_code -eq 0 ]; then
        csv_lines=$(wc -l < "$output_csv")
        echo "  ✓ Complete. Results: $csv_lines lines written to $output_csv"
    else
        echo "  ✗ Error (exit code $exit_code). Check $output_log for details."
    fi

    # Report warning summary
    if [ -f "$output_log" ] && [ -s "$output_log" ]; then
        warn_count=$(wc -l < "$output_log")
        echo "  ⚠ $warn_count warning(s) logged to $output_log"
    fi

    echo ""
done

echo "All global all-atom RMSD calculations complete."

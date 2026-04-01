#!/bin/bash

# Global TM calculation script for all proteins.
# Runs cli-global.py with --stat tm and writes per-protein CSV/log outputs.
# Usage: ./cli-global-tm.sh

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "$script_dir/.." && pwd)"

data_dir="$repo_root/data"
results_dir="$repo_root/results"
logs_dir="$results_dir/logs"

mkdir -p "$results_dir" "$logs_dir"

proteins=("PR" "IN" "RT")

echo "Starting global TM calculations for all proteins..."
echo ""

for protein in "${proteins[@]}"; do
    echo "Processing protein: $protein"

    alignment_file="$data_dir/$protein/refs/alignment.fasta"
    protein_lower=$(echo "$protein" | tr '[:upper:]' '[:lower:]')
    output_csv="$results_dir/${protein_lower}-global-tm.csv"
    output_log="$logs_dir/${protein_lower}-global-tm.stderr.log"

    if [ ! -f "$alignment_file" ]; then
        echo "  WARNING: Alignment file not found: $alignment_file"
        echo "  Skipping $protein..."
        echo ""
        continue
    fi

    echo "  Running cli-global.py..."
    python3 "$script_dir/cli-global.py" \
        --base-path "$data_dir" \
        --protein "$protein" \
        --stat tm \
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

    if [ -f "$output_log" ] && [ -s "$output_log" ]; then
        warn_count=$(wc -l < "$output_log")
        echo "  ⚠ $warn_count warning(s) logged to $output_log"
    fi

    echo ""
done

echo "All global TM calculations complete."


#!/bin/bash

# Generate positional pLDDT plots for all configured proteins.
# Expects positional pLDDT CSVs to already exist in results/.

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "$script_dir/.." && pwd)"

results_dir="$repo_root/results"
logs_dir="$results_dir/logs"
figs_dir="$results_dir/figures"

mkdir -p "$figs_dir" "$logs_dir"

proteins=("PR" "IN" "RT")

echo "Starting positional pLDDT plotting..."
echo ""

for protein in "${proteins[@]}"; do
    protein_lower=$(echo "$protein" | tr '[:upper:]' '[:lower:]')
    csv_file="$results_dir/${protein_lower}-pos-plddt.csv"
    plot_log="$logs_dir/${protein_lower}-pos-plddt.plot.stderr.log"

    echo "Processing protein: $protein"

    if [ ! -f "$csv_file" ]; then
        echo "  WARNING: CSV not found: $csv_file"
        echo "  Skipping $protein..."
        echo ""
        continue
    fi

    # Run plotting from figures dir so output PNG is written there.
    (
        cd "$figs_dir" || exit 1
        python3 "$script_dir/plot.py" \
            --csv_path "$csv_file" \
            --type plddt \
            --protein "$protein" \
            2> "$plot_log"
    )
    exit_code=$?

    if [ $exit_code -eq 0 ]; then
        echo "  ✓ Complete. Plot written to $figs_dir/${protein_lower}-pos-plddt.png"
    else
        echo "  ✗ Error (exit code $exit_code). Check $plot_log for details."
    fi

    echo ""
done

echo "All positional pLDDT plotting complete."


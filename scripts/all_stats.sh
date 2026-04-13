#!/usr/bin/env bash

script_dir="$(cd "$(dirname "$0")" && pwd)"
repo_root="$(cd "$script_dir/.." && pwd)"

results_dir="$repo_root/results"

for protein in PR IN RT; do
  python3 "$script_dir/stats.py" --csv_path "${results_dir}/${protein}-global-rmsd.csv" --type grmsd
  python3 "$script_dir/stats.py" --csv_path "${results_dir}/${protein}-global-tm.csv"   --type tm
  python3 "$script_dir/stats.py" --csv_path "${results_dir}/${protein}-pos-rmsd.csv"    --type rmsd
  python3 "$script_dir/stats.py" --csv_path "${results_dir}/${protein}-pos-plddt.csv"   --type plddt
done

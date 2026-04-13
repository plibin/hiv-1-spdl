#!/usr/bin/env bash

repo_root="$(cd "$script_dir/.." && pwd)"

results_dir="$repo_root/results"

for protein in PR IN RT; do
  python3 stats.py --csv_path "${results_dir}/${protein}-global-rmsd.csv" --type grmsd
  python3 stats.py --csv_path "${results_dir}/${protein}-global-tm.csv"   --type tm
  python3 stats.py --csv_path "${results_dir}/${protein}-pos-rmsd.csv"    --type rmsd
  python3 stats.py --csv_path "${results_dir}/${protein}-pos-plddt.csv"   --type plddt
done

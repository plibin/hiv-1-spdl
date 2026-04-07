#!/bin/bash

# Unified plotting script.
# Usage examples:
#   ./plot.sh pos rmsd
#   ./plot.sh pos plddt
#   ./plot.sh global rmsd
#   ./plot.sh global tm
#   ./plot.sh correlation

set -u

scope="${1:-pos}"
stat="${2:-rmsd}"

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "$script_dir/.." && pwd)"

results_dir="$repo_root/results"
logs_dir="$results_dir/logs"
figs_dir="$results_dir/figures"

mkdir -p "$figs_dir" "$logs_dir"

proteins=("PR" "IN" "RT")

# --- correlation scope (no stat arg needed) ---
if [ "$scope" = "correlation" ]; then
  echo "Plotting correlation (pLDDT vs RMSD) ..."
  for protein in "${proteins[@]}"; do
    protein_lower=$(echo "$protein" | tr '[:upper:]' '[:lower:]')
    rmsd_csv="$results_dir/${protein_lower}-pos-rmsd.csv"
    plddt_csv="$results_dir/${protein_lower}-pos-plddt.csv"
    plot_log="$logs_dir/${protein_lower}-correlation.plot.stderr.log"

    if [ ! -f "$rmsd_csv" ] || [ ! -f "$plddt_csv" ]; then
      echo "- $protein: missing CSV(s)"; continue
    fi

    (
      cd "$figs_dir" || exit 1
      python3 "$script_dir/plot.py" \
        --csv_path "$rmsd_csv" \
        --csv_path2 "$plddt_csv" \
        --type correlation \
        --protein "$protein" \
        --output "${protein_lower}-correlation.png" \
        2> "$plot_log"
    )
    exit_code=$?

    if [ $exit_code -eq 0 ]; then
      echo "- $protein: ok"
    else
      echo "- $protein: failed (see $plot_log)"
    fi
  done
  echo "Done."
  exit 0
fi

case "$scope" in
  pos)
    case "$stat" in
      rmsd)
        csv_suffix="pos-rmsd"
        plot_type="rmsd"
        ;;
      plddt)
        csv_suffix="pos-plddt"
        plot_type="plddt"
        ;;
      *)
        echo "Invalid stat for scope 'pos': $stat (use: rmsd|plddt)"
        exit 2
        ;;
    esac
    ;;
  global)
    case "$stat" in
      rmsd)
        csv_suffix="global-rmsd"
        plot_type="grmsd"
        ;;
      tm)
        csv_suffix="global-tm"
        plot_type="tm"
        ;;
      *)
        echo "Invalid stat for scope 'global': $stat (use: rmsd|tm)"
        exit 2
        ;;
    esac
    ;;
  *)
    echo "Invalid scope: $scope (use: pos|global)"
    exit 2
    ;;
esac

echo "Plotting $scope/$stat ..."

for protein in "${proteins[@]}"; do
  protein_lower=$(echo "$protein" | tr '[:upper:]' '[:lower:]')
  csv_file="$results_dir/${protein_lower}-${csv_suffix}.csv"
  plot_log="$logs_dir/${protein_lower}-${csv_suffix}.plot.stderr.log"

  [ -f "$csv_file" ] || { echo "- $protein: missing CSV"; continue; }

  (
    cd "$figs_dir" || exit 1
    python3 "$script_dir/plot.py" \
      --csv_path "$csv_file" \
      --type "$plot_type" \
      --protein "$protein" \
      2> "$plot_log"
  )
  exit_code=$?

  if [ $exit_code -eq 0 ]; then
    echo "- $protein: ok"
  else
    echo "- $protein: failed (see $plot_log)"
  fi
done

echo "Done."


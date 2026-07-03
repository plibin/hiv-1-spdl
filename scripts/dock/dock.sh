#!/usr/bin/env bash
set -euo pipefail

#e.g.: ./dock.sh ligand.sdf clean_pr.pdb 22 A:25:OD1,OD2 B:25:OD1,OD2

# mk_prepare_* comes from the meeko package (pip install meeko)
# shift x removes the first x command-line arguments from $@

#conda install -c conda-forge vina

ligand="$1"; protein="$2"; size="${3}"
shift 3

ligbase="$(basename "$ligand")"; ligbase="${ligbase%.*}"
protbase="$(basename "$protein")"; protbase="${protbase%.*}"

out_dir="$(dirname "$protein")/results_${protbase}_${ligbase}"
prep_dir="${out_dir}/prep"
dock_dir="${out_dir}/docking"

mkdir -p "$prep_dir" "$dock_dir"

receptor_pdbqt="${prep_dir}/${protbase}.receptor.pdbqt"
ligand_pdbqt="${prep_dir}/${ligbase}.ligand.pdbqt"

read cx cy cz < <(python pocket_center.py "$protein" "$@")

protonated_pdb="${prep_dir}/${protbase}.protonated.pdb"
# ASSUMPTION: Default pH 7.4 protonation deprotonates both Asp25 residues. Unrealistic because HIV PR mechanism requires an asymmetric dyad.
# Fix: Use PROPKA to calculate local pKa shifts and assign correct protonation, and reconstruct missing heavy atoms.
# Reference: https://pmc.ncbi.nlm.nih.gov/articles/PMC2660780/
pdb2pqr --ff=AMBER --titration-state-method=propka --with-ph=7.4 --pdb-output "$protonated_pdb" "$protein" "${prep_dir}/${protbase}.pqr"
mv *.propka "$prep_dir/" 2>/dev/null || true

mk_prepare_receptor.py -i "$protonated_pdb" -o "${prep_dir}/${protbase}.receptor" -p -a > "${prep_dir}/${protbase}.receptor.log"

mv "${prep_dir}/${protbase}.receptor.pdbqt" "$receptor_pdbqt"

mk_prepare_ligand.py -i "$ligand" -o "$ligand_pdbqt" > "$ligand_pdbqt.log"

csv_file="${out_dir}/${protbase}_${ligbase}.scores.csv"

for i in {1..2}; do
  seed=$i
  current_out="${dock_dir}/${protbase}_${ligbase}_seed${seed}.docked.pdbqt"
  current_log="${dock_dir}/${protbase}_${ligbase}_seed${seed}.vina.log"
  
  vina --receptor "$receptor_pdbqt" --ligand "$ligand_pdbqt" \
    --center_x "$cx" --center_y "$cy" --center_z "$cz" \
    --size_x "$size" --size_y "$size" --size_z "$size" \
    --seed "$seed" \
    --out "$current_out" > "$current_log" 2>&1
    
  top_score=$(grep "^REMARK VINA RESULT" "$current_out" | head -n 1 | awk '{print $4}')
  echo "${protbase},${seed},${top_score}" 
done

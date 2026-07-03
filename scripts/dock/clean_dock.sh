f=$1

if [ -f "$f" ]; then
  python clean_pdb.py "$f" "$f.clean"
  ./dock.sh ./Darunavir.sdf "$f.clean" 22 A:25:OD1,OD2 B:25:OD1,OD2
fi

echo 'name,seed,score' > AF2.docking.csv
echo 'name,seed,score' > AF3.docking.csv
echo 'name,seed,score' > refs.docking.csv

for f in ./dock_data/refs/*.pdb; do
  name="$(basename "$f")"
  ./clean_dock.sh "./dock_data/AlphaFold2/$name" >> AF2.docking.csv 
  ./clean_dock.sh "./dock_data/AlphaFold3/$name" >> AF3.docking.csv 
  ./clean_dock.sh "./dock_data/refs/$name" >> refs.docking.csv 
done

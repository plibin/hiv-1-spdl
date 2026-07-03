cp -r ../../data .

data=./data

#in the refs, we have experimental PDBs, so there might be quality issues,
#so we use these to get our initial selection, we will use this selection to run the dockings for the algorithms
mkdir -p dock_data/refs
find $data/PR/refs/ -name "*.pdb" -type f -exec python select_pdb.py {} ./dock_data/refs/ \;

mkdir -p dock_data/AlphaFold2
find $data/PR/AlphaFold2/ -name "*.pdb" -type f -exec python select_pdb.py {} ./dock_data/AlphaFold2/ \;

mkdir -p dock_data/AlphaFold3
find $data/PR/AlphaFold3/ -name "*.pdb" -type f -exec python select_pdb.py {} ./dock_data/AlphaFold3/ \;

mkdir -p dock_data/ESMFold
find $data/PR/ESMFold/ -name "*.pdb" -type f -exec python select_pdb.py {} ./dock_data/ESMFold/ \;

cd dock_data/ESMFold
for f in *.pdb; do mv "$f" "$(echo "${f:0:4}" | tr A-Z a-z).pdb"; done

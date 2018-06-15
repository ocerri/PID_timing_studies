unset LD_LIBRARY_PATH
unset PYTHONHOME
unset PYTHONPATH
source /cvmfs/sft.cern.ch/lcg/releases/LCG_88/gcc/4.9.3/x86_64-slc6/setup.sh
source /cvmfs/sft.cern.ch/lcg/releases/LCG_88/Python/2.7.13/x86_64-slc6-gcc49-opt/Python-env.sh
source /cvmfs/sft.cern.ch/lcg/releases/LCG_88/ROOT/6.08.06/x86_64-slc6-gcc49-opt/bin/thisroot.sh
source /cvmfs/sft.cern.ch/lcg/releases/LCG_88/scipy/0.18.1/x86_64-slc6-gcc49-opt/scipy-env.sh

mkdir job
cd job

cp /afs/cern.ch/work/s/selvaggi/private/Delphes/libDelphes.so .

INPUT_FILE=${1}
OUTPUT_FILE=${2}

cp ${INPUT_FILE} in.root

python /afs/cern.ch/work/s/selvaggi/private/Delphes/examples/produceTree_MLWS.py in.root out.root
cp out.root ${OUTPUT_FILE}


#!/bin/bash

source /afs/cern.ch/user/o/ocerri/cernbox/PID_timing_studies/bin/setup_lxp.sh

DELPHESCARD=${1}
PYTHIACARD=${2}
PTHAT_min=${3}
PTHAT_max=${4}
NEVENTS=${5}
SEED=${6}
OUTPUTFILE=${7}
STOPMASS=${8}

mkdir job_${SEED}
cd job_${SEED}

cp $main_dir/_hepmc/MinBias_pp14TeV_100k.pileup ./
cp $main_dir/Generation/delphes/cards/* ./

cp $main_dir/Generation/pythia/cards/* ./
sed -i -e '37s/.*/   1000006     '$STOPMASS'    # ~t_1/' split_stop.slha
mkdir cards
mv split_stop.slha cards/split_stop.slha

echo 'Random:setSeed = on' >> ${PYTHIACARD}
echo 'Random:seed = '${SEED} >> ${PYTHIACARD}
echo 'PhaseSpace:pTHatMin = '${PTHAT_min} >> ${PYTHIACARD}
echo 'PhaseSpace:pTHatMax = '${PTHAT_max} >> ${PYTHIACARD}
echo 'Main:numberOfEvents = '${NEVENTS} >> ${PYTHIACARD}

$Delphes/DelphesPythia8 ${DELPHESCARD} ${PYTHIACARD} ${OUTPUTFILE}.root &> ${OUTPUTFILE}.log
# cp out.root ${OUTPUTFILE}.root

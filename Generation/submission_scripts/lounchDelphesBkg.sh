
#!/bin/bash

source /afs/cern.ch/user/o/ocerri/cernbox/PID_timing_studies/bin/setup_lxp7.sh

DELPHESCARD=CMS_PhaseII_trackNtime_lightTree.tcl
PYTHIACARD=pp2HardQCD.pythia
PTHAT_min=200
PTHAT_max=600
NEVENTS=40000
SEED=${1}
OUTPUTFILE=/afs/cern.ch/user/o/ocerri/cernbox/PID_timing_studies/_root/jobs/pp2HardQCD_PU100_pTHat${PTHAT_min}-${PTHAT_max}_${SEED}

cd $Delphes

mkdir job_${SEED}

DELPHESCARD=$main_dir/Generation/delphes/cards/${DELPHESCARD}
cp $main_dir/Generation/pythia/cards/* ./job_${SEED}/
PYTHIACARD=job_${SEED}/${PYTHIACARD}


echo 'Random:setSeed = on' >> ${PYTHIACARD}
echo 'Random:seed = '${SEED} >> ${PYTHIACARD}
echo 'PhaseSpace:pTHatMin = '${PTHAT_min} >> ${PYTHIACARD}
echo 'PhaseSpace:pTHatMax = '${PTHAT_max} >> ${PYTHIACARD}
echo 'Main:numberOfEvents = '${NEVENTS} >> ${PYTHIACARD}

./DelphesPythia8 ${DELPHESCARD} ${PYTHIACARD} ${OUTPUTFILE}.root >> ${OUTPUTFILE}.log

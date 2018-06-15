#!/bin/bash

source /cvmfs/sft.cern.ch/lcg/releases/LCG_88/gcc/4.9.3/x86_64-slc6/setup.sh
source /cvmfs/sft.cern.ch/lcg/releases/LCG_88/Python/2.7.13/x86_64-slc6-gcc49-opt/Python-env.sh
source /cvmfs/sft.cern.ch/lcg/releases/LCG_88/ROOT/6.08.06/x86_64-slc6-gcc49-opt/bin/thisroot.sh

source /cvmfs/sft.cern.ch/lcg/releases/LCG_88/Python/2.7.13/x86_64-slc6-gcc49-opt/Python-env.sh
source /cvmfs/sft.cern.ch/lcg/releases/LCG_88/ROOT/6.08.06/x86_64-slc6-gcc49-opt/ROOT-env.sh
source /cvmfs/sft.cern.ch/lcg/releases/LCG_88/scipy/0.18.1/x86_64-slc6-gcc49-opt/scipy-env.sh


#Getting pythia8
source /afs/cern.ch/user/o/ocerri/work/programs_180316/pythia8230/user_setup.sh
export Pythia=$PYTHIA8

#Getting Delphes431
export Delphes=/afs/cern.ch/user/o/ocerri/work/programs_180316/Delphes
cd $Delphes
. DelphesEnv.sh
cd -

#Adding the to the python path the lib folder
export PYTHONPATH=/afs/cern.ch/user/o/ocerri/cernbox/PID_timing_studies/lib:$PYTHONPATH

export main_dir=/afs/cern.ch/user/o/ocerri/cernbox/PID_timing_studies
alias cd_main="cd $main_dir"

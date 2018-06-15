#!/bin/bash

#Getting root 6.10.08
. /cvmfs/sft.cern.ch/lcg/contrib/gcc/4.8/x86_64-centos7-gcc48-opt/setup.sh
. /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.10.08/x86_64-centos7-gcc48-opt/root/bin/thisroot.sh

#Getting pythia8
source /afs/cern.ch/user/o/ocerri/work/programs_171206/pythia8230/user_setup.sh
export Pythia=$PYTHIA8

#Getting Delphes431
export Delphes=/afs/cern.ch/user/o/ocerri/work/programs_171206/Delphes-3.4.1
cd $Delphes
. DelphesEnv.sh
cd -

#Adding the to the python path the lib folder
export PYTHONPATH=/afs/cern.ch/user/o/ocerri/cernbox/PID_timing_studies/lib:$PYTHONPATH

export main_dir=/afs/cern.ch/user/o/ocerri/cernbox/PID_timing_studies
alias cd_main="cd $main_dir"

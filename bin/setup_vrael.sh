#!/bin/bash

#Getting pythia8
source $HOME/programs/pythia8230/user_setup.sh
export Pythia=$PYTHIA8

#Getting Delphes431
export Delphes=$HOME/programs/Delphes-3.4.1
cd $Delphes
. DelphesEnv.sh
cd -

#Adding the to the python path the lib folder
export PYTHONPATH=$HOME/cernbox/PID_timing_studies/lib:$PYTHONPATH
export PATH=$HOME/cernbox/PID_timing_studies/lib:$PATH

export main_dir=$HOME/cernbox/PID_timing_studies
alias cd_main="cd $main_dir"

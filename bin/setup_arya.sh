#!/bin/bash

#Getting Delphes431
export Delphes=$HOME/programs/Delphes-3.4.1
cd $Delphes
. DelphesEnv.sh
cd -

#Adding the to the python path the lib folder
export PYTHONPATH=$HOME/cernbox/PID_timing_studies/lib:$PYTHONPATH

export main_dir=$HOME/cernbox/PID_timing_studies
alias cd_main="cd $main_dir"

#!/bin/bash

python src/extract_flat_evts_tree.py /Users/olmo/cernbox/PID_timing_studies/_root/jobs_PU140_poiss/pp2HardQCD_PU140_pTHat150-185/*.root
python src/extract_flat_evts_tree.py /Users/olmo/cernbox/PID_timing_studies/_root/jobs_PU140_poiss/pp2HardQCD_PU140_pTHat300-600/*.root &> log300.log &

echo $Delphes
cd $Delphes

make -j8 DelphesHepMC
rm $main_dir/Generation/delphes/log.log
rm -rf ~/Desktop/debug/*
rm $main_dir/_root/$1.root

./DelphesHepMC $main_dir/Generation/delphes/cards/CMS_PhaseII_trackNtime.tcl $main_dir/_root/$1.root $main_dir/_hepmc/gg2stopantistop_RHad_10k_v1.hepmc >> $main_dir/Generation/delphes/log.log

cd -

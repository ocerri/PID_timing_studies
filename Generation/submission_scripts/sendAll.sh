rm -rf /eos/user/s/selvaggi/MLWS/WW_Pt5000_0PU; python submitDelphesJobs.py --pt 5000. --pycard configBoostedw.cmnd --dcard FCC/FCChh.tcl --outdir /eos/user/s/selvaggi/MLWS/WW_Pt5000_0PU --njobs 100 --nev 1500 --queue 1nd --memory 16000. --disk 8000.
rm -rf /eos/user/s/selvaggi/MLWS/QCD_Pt5000_0PU; python submitDelphesJobs.py --pt 5000. --pycard configBoostedqcd_weak.cmnd --dcard FCC/FCChh.tcl --outdir /eos/user/s/selvaggi/MLWS/QCD_Pt5000_0PU --njobs 100 --nev 1500 --queue 1nd --memory 16000. --disk 8000.

python submitAnalysisJobs.py -i /eos/user/s/selvaggi/MLWS/WW_Pt5000_0PU -o /eos/user/s/selvaggi/MLWS/WW_Pt5000_0PU/w
python submitAnalysisJobs.py -i /eos/user/s/selvaggi/MLWS/QCD_Pt5000_0PU -o /eos/user/s/selvaggi/MLWS/QCD_Pt5000_0PU/qcd


rm -rf /afs/cern.ch/work/s/selvaggi/public/4MLWS/qcd/* /afs/cern.ch/work/s/selvaggi/public/4MLWS/w/*

cp -r /eos/user/s/selvaggi/MLWS/WW_Pt5000_0PU/w/out/* /afs/cern.ch/work/s/selvaggi/public/4MLWS/w/
cp -r /eos/user/s/selvaggi/MLWS/QCD_Pt5000_0PU/qcd/out/* /afs/cern.ch/work/s/selvaggi/public/4MLWS/qcd



rm -rf /eos/user/s/selvaggi/MLWS/QCD_Pt5000_0PU; python submitDelphesJobs.py --pt 5000. --pycard configBoostedqcd_weak.cmnd --dcard FCC/FCChh_MLWS.tcl --outdir /eos/user/s/selvaggi/MLWS/QCD_Pt5000_0PU --njobs 100 --nev 1000 --queue 1nd --memory 16000. --disk 8000.
rm -rf /eos/user/s/selvaggi/MLWS/WW_Pt5000_0PU; python submitDelphesJobs.py --pt 5000. --pycard configBoostedw.cmnd --dcard FCC/FCChh_MLWS.tcl --outdir /eos/user/s/selvaggi/MLWS/WW_Pt5000_0PU --njobs 100 --nev 1000 --queue 1nd --memory 16000. --disk 8000.


python submitAnalysisJobs.py -i /eos/user/s/selvaggi/MLWS/QCD_Pt5000_0PU -o /eos/user/s/selvaggi/MLWS/QCD_Pt5000_0PU/qcd
python submitAnalysisJobs.py -i /eos/user/s/selvaggi/MLWS/QCD_Pt5000_0PU -o /eos/user/s/selvaggi/MLWS/QCD_Pt5000_0PU/qcd



cp -r /eos/user/s/selvaggi/MLWS/QCD_Pt5000_0PU/qcd/out/* /afs/cern.ch/work/s/selvaggi/public/4MLWS/qcd

import numpy as np
from root_numpy import array2root, root2array

file_path = "/Users/olmo/cernbox/PID_timing_studies/_root/jobs_PU140_poiss/pp2HardQCD_PU140_pTHat200-INF/pp2HardQCD_PU140_pTHat200-INF_950k_tks_flat.root"
file_out = "/Users/olmo/cernbox/PID_timing_studies/_root/histos/RecoMassTree_bkg.root"

arr = root2array(file_path, treename="T", branches="M_reco", selection="M_reco>30 && pt>100 && vtx_SumPT>370")
arr = np.atleast_2d(arr.flatten()).T
arr = np.array(arr, dtype=[('Mass', np.float32)])

array2root(arr, file_out, treename='T', mode='RECREATE')

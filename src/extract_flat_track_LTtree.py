import numpy as np
import ROOT as rt
from root_numpy import array2root
import os, sys

c_light = 2.99792458E8 #m/s

current_dir = os.getcwd()
os.chdir(os.environ['Delphes'])
rt.gSystem.Load("libDelphes");
rt.gInterpreter.Declare('#include "classes/DelphesClasses.h"');
rt.gInterpreter.Declare('#include "external/ExRootAnalysis/ExRootTreeReader.h"');


file_path = list(sys.argv[1:])
print len(sys.argv)-1

chain = rt.TChain('Delphes')
for path in file_path:
    chain.Add(path)

treeReader = rt.ExRootTreeReader(chain)
numberOfEntries = treeReader.GetEntries()
print numberOfEntries

branches_names = ['Event', 'Track','Vertex4D']
branch = {}

for b in branches_names:
    branch[b] = treeReader.UseBranch(b)



tks_names = ['Nev','tof_reco', 'PID', 'tof_gen',
             'P_reco', 'vtx_SumPT2', 'vtx_NDOF','vtx_SumPT', 'i_vtx',
             'Zout','Tout','pt','ctgtheta','phi','d0','dz', 'L',
             'sigma_pt', 'sigma_d0', 'sigma_dz', 'sigma_Tin',
             'M_reco', 'beta_reco'
            ]
tks_att = []


for i in range(numberOfEntries):
    if i%50000==0:
        print i
    treeReader.ReadEntry(i)

    Nev = i
    for j in range(branch['Track'].GetEntries()):
        track = branch['Track'].At(j)
        if track.VertexIndex == -1: continue

        i_vtx = track.VertexIndex
        vtx = branch['Vertex4D'].At(i_vtx)

        tof_gen = 1E12*track.tof_gen

        tof_reco = 1E12*track.tof_reco


        PID = track.PID

        P_reco = track.P

        beta_reco_fromL = track.L*1E-3/(c_light*1E-12*tof_reco)

        M_reco = track.Mass

        tks_att.append((Nev,tof_reco, PID, tof_gen,
                        P_reco,
                        vtx.SumPT2, vtx.NDF, vtx.SumPT, i_vtx,
                        track.ZOuter, 1E12*track.TOuter, track.PT, track.CtgTheta, track.Phi,
                        track.D0, track.DZ, track.L,
                        track.ErrorPT, track.ErrorD0, track.ErrorDZ, 1E12*vtx.ErrorT,
                        M_reco, beta_reco_fromL
                        ))



tks = np.array(tks_att, dtype=zip(tks_names,['<f8']*len(tks_names)))
print tks.size
#
array2root(tks, file_path[0].replace('1.root', '{}k_tks_flat.root'.format(numberOfEntries/1000)), treename='T', mode='RECREATE')

os.chdir(current_dir)

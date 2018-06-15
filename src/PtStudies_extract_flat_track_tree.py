import numpy as np
import ROOT as rt
from root_numpy import array2root
import os, sys

c_light = 2.99792458E8 #m/s

current_dir = os.getcwd()
os.chdir('/Users/olmo/programs/Delphes-3.4.1')
rt.gSystem.Load("libDelphes");
rt.gInterpreter.Declare('#include "classes/DelphesClasses.h"');
rt.gInterpreter.Declare('#include "external/ExRootAnalysis/ExRootTreeReader.h"');


file_path = list(sys.argv[1:])
print sys.argv
print file_path

chain = rt.TChain('Delphes')
for path in file_path:
    chain.Add(path)

treeReader = rt.ExRootTreeReader(chain)
numberOfEntries = treeReader.GetEntries()

branches_names = ['Event', 'Particle', 'Track','Vertex4D', 'GenVertex']
branch = {}

for b in branches_names:
    branch[b] = treeReader.UseBranch(b)



# tks_names = ['Nev','DzOF','tof_reco', 'IsPU', 'PID', 'M_gen', 'tof_gen', 'P_reco', 'P_gen', 'CtgTheta', 'i_vtx',
#              'vtx_SumPT2', 'vtx_NDOF','vtx_SumPT',
#              'Xout','Yout','Zout','Tout','pt','ctgtheta','phi','d0','dz', 'L',
#              'sigma_Tout', 'sigma_pt', 'sigma_d0', 'sigma_dz', 'sigma_Zout', 'sigma_Tin',
#              'M_reco', 'beta_reco', 'beta_MC', 'error_Tin', 'error_Zin', 'Tin_MC', 'Zin_MC',
#              'Td'
#             ]
tks_names = ['Nev', 'PID', 'M_gen', 'pt_gen', 'eta_gen', 'phi_gen', 'beta_gen',
             'M_reco', 'pt_reco', 'eta_reco', 'phi_reco', 'beta_reco',
             'vtx_idx', 'vtx_SumPT'
            ]
tks_att = []


for i in range(numberOfEntries):
    treeReader.ReadEntry(i)

    Nev = i

    for j in range(branch['Track'].GetEntries()):
        track = branch['Track'].At(j)
        p = track.Particle.GetObject()

        if abs(p.PID) < 100000 : continue

        beta_gen = p.P/p.E
        beta_reco = track.L*1E-3/(c_light * track.tof_reco)

        i_vtx = track.VertexIndex
        vtx_SumPT = -1
        if track.VertexIndex >= 0:
            vtx_SumPT = branch['Vertex4D'].At(i_vtx).SumPT

        tks_att.append((Nev, p.PID, p.Mass, p.PT, p.Eta, p.Phi, beta_gen,
                        track.Mass, track.PT, track.Eta, track.Phi, beta_reco,
                        i_vtx, vtx_SumPT
                        ))



tks = np.array(tks_att, dtype=zip(tks_names,['<f8']*len(tks_names)))
print tks.size

array2root(tks, file_path[0].replace('.root', '_BSMtks_flat.root'), treename='T', mode='RECREATE')

os.chdir(current_dir)

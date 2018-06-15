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



tks_names = ['Nev','DzOF','tof_reco', 'IsPU', 'PID', 'M_gen', 'tof_gen', 'P_reco', 'P_gen', 'CtgTheta', 'i_vtx',
             'vtx_SumPT2', 'vtx_NDOF','vtx_SumPT',
             'Xout','Yout','Zout','Tout','pt','ctgtheta','phi','d0','dz', 'L',
             'sigma_Tout', 'sigma_pt', 'sigma_d0', 'sigma_dz', 'sigma_Zout', 'sigma_Tin',
             'M_reco', 'beta_reco', 'beta_MC', 'error_Tin', 'error_Zin', 'Tin_MC', 'Zin_MC',
             'Td'
            ]
tks_att = []


for i in range(numberOfEntries):
    treeReader.ReadEntry(i)

    Nev = i

    for j in range(branch['Track'].GetEntries()):
        track = branch['Track'].At(j)
        if track.VertexIndex == -1: continue

        i_vtx = track.VertexIndex
        vtx = branch['Vertex4D'].At(i_vtx)

        p = track.Particle.GetObject()
        beta_MC = p.P/p.E

        tof_gen = 1E12*track.tof_gen

        tof_reco = 1E12*track.tof_reco

        error_Tin = 1E12*(vtx.T - p.T)
        error_Zin = track.Z - p.Z

        Zin_MC = p.Z
        Tin_MC = 1E12*p.T

        DzOF = track.ZOuter - vtx.Z

        IsPU = p.IsPU

        M = p.Mass

        PID = p.PID

        P_reco = track.P
        P_gen = p.P

        CtgTheta = track.CtgTheta

        beta_reco_fromZ = np.sqrt(1 + 1./np.square(CtgTheta) ) * np.abs(DzOF)*1E-3/(c_light*1E-12*tof_reco)
        beta_reco_fromL = track.L*1E-3/(c_light*1E-12*tof_reco)

        # M_reco = P*np.sqrt(1./np.square(np.minimum(beta_reco,np.ones_like(beta_reco))) - 1)
        M_reco = track.Mass

        tks_att.append((Nev,DzOF,tof_reco, IsPU, PID, M, tof_gen, P_reco, P_gen, CtgTheta, i_vtx,
                        vtx.SumPT2, vtx.NDF, vtx.SumPT,
                        track.XOuter, track.YOuter, track.ZOuter, 1E12*track.TOuter, track.PT, track.CtgTheta, track.Phi, track.D0, track.DZ,
                        track.L,
                        1E12*track.ErrorTOut, track.ErrorPT, track.ErrorD0, track.ErrorDZ, track.ErrorZOut, 1E12*vtx.ErrorT,
                        M_reco, beta_reco_fromL, beta_MC, error_Tin, error_Zin, Tin_MC, Zin_MC, 1E12*track.Td
                        ))



tks = np.array(tks_att, dtype=zip(tks_names,['<f8']*len(tks_names)))
print tks.size

array2root(tks, file_path[0].replace('.root', '_tks_flat.root'), treename='T', mode='RECREATE')

os.chdir(current_dir)

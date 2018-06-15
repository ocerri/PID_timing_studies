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


# gvtx_names = ['Nev','z','t', 'IsPU']
# gvtx_att = []

rvtx_names = ['Nev','z','t', 'dz', 'dt', 'ndf', 'sumpt2', 'sum_pt', 'N_Rhad', 'max_pt_tk']
rvtx_att_global = []


for i in range(numberOfEntries):
    treeReader.ReadEntry(i)

    # if i>2:
    #     break

    # for j in range(branch['GenVertex'].GetEntries()):
    #     vtx = branch['GenVertex'].At(j)
    #
    #     p = vtx.Constituents.At(0)
    #     ispu = p.IsPU
    #
    #     gvtx_att.append((i, vtx.Z, vtx.T*1E12, ispu))

    rvtx_att = []
    for j in range(branch['Vertex4D'].GetEntries()):
        vtx = branch['Vertex4D'].At(j)

        rvtx_att.append([i, vtx.Z, vtx.T*1E12, vtx.ErrorZ, vtx.ErrorT*1E12, vtx.NDF, vtx.SumPT2, vtx.SumPT, 0, 0])

    # print rvtx_att

    for j in range(branch['Track'].GetEntries()):
        track = branch['Track'].At(j)
        if track.VertexIndex >= 0:
            k = track.VertexIndex
            pt = track.PT

            if pt > rvtx_att[k][9]:
                rvtx_att[k][9] = pt

            if abs(track.PID)>5000:
                rvtx_att[k][8] += 1
    # print rvtx_att
    # print rvtx_att[0]

    for j in range(len(rvtx_att)):
        rvtx_att_global.append(tuple(rvtx_att[j]))

    # print rvtx_att_global


# gen_vtx = np.array(gvtx_att, dtype=zip(gvtx_names,['<f8']*len(gvtx_names)))
reco_vtx = np.array(rvtx_att_global, dtype=zip(rvtx_names,['<f8']*len(rvtx_names)))

print reco_vtx.shape

array2root(reco_vtx, file_path[0].replace('.root', '_vtx_flat.root'), treename='T', mode='RECREATE')

os.chdir(current_dir)

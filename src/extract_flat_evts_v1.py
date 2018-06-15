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
    if '_tks_' in path:
        continue
    chain.Add(path)

treeReader = rt.ExRootTreeReader(chain)
numberOfEntries = treeReader.GetEntries()

branches_names = ['Event', 'Track','Vertex4D']
branch = {}

for b in branches_names:
    branch[b] = treeReader.UseBranch(b)



features_names = ['Nev', 'mh', 'pth', 'PIDh', 'ml', 'ptl', 'PIDl', 'vtx_SumPT'
            ]
evt_features = []


for i in range(numberOfEntries):
    if i%50000 == 0:
        print i
    treeReader.ReadEntry(i)
    N_tracks = branch['Track'].GetEntries()
    if N_tracks < 1: continue

    m_h = 0
    i_h = 0
    for j in range(N_tracks):
        mass = branch['Track'].At(j).Mass
        if mass > m_h:
            m_h = mass
            i_h = j

    t_h = branch['Track'].At(i_h)
    i_vtx = t_h.VertexIndex
    vtx = branch['Vertex4D'].At(i_vtx)

    m_l = -1
    i_l = -1
    for j in range(N_tracks):
        if j == i_h or branch['Track'].At(j).VertexIndex != i_vtx:
            continue

        mass = branch['Track'].At(j).Mass
        if mass > m_l:
            m_l = mass
            i_l = j

    if m_l == -1:
        evt_features.append((i, t_h.Mass, t_h.PT, t_h.PID, -1, -1, 0, vtx.SumPT))
    else:
        t_l = branch['Track'].At(i_l)
        evt_features.append((i, t_h.Mass, t_h.PT, t_h.PID, t_l.Mass, t_l.PT, t_l.PID, vtx.SumPT))



tks = np.array(evt_features, dtype=zip(features_names,['<f8']*len(features_names)))
print tks.size

array2root(tks, file_path[0].replace('1.root', '{}k_evt_flat.root'.format(numberOfEntries/1000)), treename='T', mode='RECREATE')


os.chdir(current_dir)

import numpy as np
import ROOT as rt
from root_numpy import array2root
import os, sys, re

c_light = 2.99792458E8 #m/s

def GetXsec(path):
    xsec_files = os.path.dirname(path) + '/xsec.txt'
    fxsec = open(xsec_files)
    err = []
    xsec = []
    for l in fxsec.readlines():
        out = re.search(r'[0-9]+.[0-9]+e-[0-9]+  [0-9]+.[0-9]+e-[0-9]+', l)
        out =  str(out.group(0))
        out = out.split('  ')
        xsec.append(float(out[0])*1e12)
        err.append(float(out[1])*1e12)

    err2 = np.square(np.array(err))
    xsec = np.array(xsec)
    xsec = np.sum(xsec*err2/np.sum(err2))
    err = 1./np.sqrt(np.sum(1./err2))

    print 'XSec = {:.1e} +/- {:.1e} fb'.format(xsec, err)

    return (xsec, err)


current_dir = os.getcwd()
os.chdir('/Users/olmo/programs/Delphes-3.4.1')
rt.gSystem.Load("libDelphes");
rt.gInterpreter.Declare('#include "classes/DelphesClasses.h"');
rt.gInterpreter.Declare('#include "external/ExRootAnalysis/ExRootTreeReader.h"');


file_path = list(sys.argv[1:])

for p in file_path:
    if '_evt' in p:
        file_path.remove(p)

print 'From dir', os.path.dirname(file_path[0])
print [os.path.basename(path) for path in file_path]

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



features_names = ['Nev', 'mh', 'pth', 'PIDh', 'ml', 'ptl', 'PIDl', 'vtx_SumPT','max_vtx_SumPT', 'same_vtx_idx', 'exist_RHad_trk'
            ]
evt_features = []

print 'Reading', numberOfEntries, 'events'
for i in range(numberOfEntries):
    if i%50000 == 0:
        print i
    treeReader.ReadEntry(i)

    N_vtx = branch['Vertex4D'].GetEntries()
    max_vtx_SumPT = 0
    idx_max_SumPT = 0
    for j in range(N_vtx):
        vtx = branch['Vertex4D'].At(j)
        if vtx.SumPT > max_vtx_SumPT:
            max_vtx_SumPT = vtx.SumPT
            idx_max_SumPT = j

    N_tracks = branch['Track'].GetEntries()
    if N_tracks < 1:
        evt_features.append((i, -1, 0, 0, -1, -1, 0, -1, max_vtx_SumPT, -1, 0))
        continue


    m_h = 0
    i_h = 0
    exist_RHad_trk = 0
    for j in range(N_tracks):
        PID = branch['Track'].At(j).PID
        if abs(PID) > 1e5:
            exist_RHad_trk += 1

        mass = branch['Track'].At(j).Mass
        if mass > m_h:
            m_h = mass
            i_h = j

    t_h = branch['Track'].At(i_h)
    i_vtx = t_h.VertexIndex

    same_vtx_idx = (i_vtx == idx_max_SumPT)

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
        evt_features.append((i, t_h.Mass, t_h.PT, t_h.PID, -1, -1, 0, vtx.SumPT, max_vtx_SumPT, same_vtx_idx, exist_RHad_trk))
    else:
        t_l = branch['Track'].At(i_l)
        evt_features.append((i, t_h.Mass, t_h.PT, t_h.PID, t_l.Mass, t_l.PT, t_l.PID, vtx.SumPT, max_vtx_SumPT, same_vtx_idx, exist_RHad_trk))



evts = np.array(evt_features, dtype=zip(features_names,['<f8']*len(features_names)))
print 'Events saved: ',evts.size

name_out = file_path[0].replace('1.root', 'evtV2_{}k.root'.format(numberOfEntries/1000))
print 'Saving in file:', name_out
array2root(evts, name_out, treename='T', mode='RECREATE')

f = rt.TFile(name_out, 'UPDATE')

aux = rt.TVectorD(2)
# aux.SetName('XSec')
# aux.SetTitle('XSec')
xsec = GetXsec(name_out)
aux[0] = xsec[0]
aux[1] = xsec[1]
aux.Write('XSec')

aux = rt.TVectorD(1)
# aux.SetName('NEvts')
# aux.SetTitle('NEvts')
aux[0] = numberOfEntries
aux.Write('NEvts')

f.Close()


os.chdir(current_dir)

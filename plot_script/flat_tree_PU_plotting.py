import numpy as np
import ROOT as rt
import root_numpy as rtnp
from glob import glob
from histo_utilities import create_TH1D, create_TH2D
from cebefo_style import cebefo_style
import re

cebefo_style()

file_list = glob('_root/gg2Rhad_PU*_tks_flat.root')
color = [1,2,4,8]
c_out = rt.TCanvas('c_Tin', 'c_Tin', 800, 600)
histo_list = []
leg = rt.TLegend(0.2,0.3)
for i, fname in enumerate(file_list):
    print fname
    res = re.search('PU[0-9]+',fname)
    PU = int(res.group(0)[2:])
    particles = rtnp.root2array(fname)

    h = create_TH1D(particles['error_Tin'], name='ht{}'.format(i),
                    title = 'Inner time resolution for different PU', binning = [151, -100,100])
    h.SetMarkerColor(color[i])
    h.SetMarkerStyle(24)
    h.SetMarkerSize(0.5)
    h.SetLineColor(color[i])
    h.Scale(1./h.Integral())
    h.SetStats(0)
    h.SetXTitle('T_{in}^{MC}-T_{in}^{reco} [ps]')
    h.Draw('sames')
    leg.AddEntry(h.GetName(), 'PU {}, std = {:.1f}'.format(PU, h.GetRMS()), 'lep')

    histo_list.append(h)

c_out.SetLogy()
leg.Draw()
c_out.Update()


c_Zin = rt.TCanvas('c_Zin', 'c_Zin', 800, 600)
leg_Z = rt.TLegend(0.2,0.3)
for i, fname in enumerate(file_list):
    print fname
    res = re.search('PU[0-9]+',fname)
    PU = int(res.group(0)[2:])
    particles = rtnp.root2array(fname)

    h = create_TH1D(particles['error_Zin'], name='hz{}'.format(i),
                    title = 'Production Z resolution for different PU', binning = [151, -1,1])
    h.SetMarkerColor(color[i])
    h.SetMarkerStyle(24)
    h.SetMarkerSize(0.5)
    h.SetLineColor(color[i])
    h.Scale(1./h.Integral())
    h.SetStats(0)
    h.SetXTitle('Z_{in}^{MC}-Z_{in}^{reco} [mm]')
    h.Draw('sames')
    leg_Z.AddEntry(h.GetName(), 'PU {}, std = {:.3f}'.format(PU, h.GetRMS()), 'lep')

    histo_list.append(h)

c_Zin.SetLogy()
leg_Z.Draw()
c_Zin.Update()

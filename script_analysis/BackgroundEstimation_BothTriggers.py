import numpy as np
import ROOT as rt
# from root_numpy import root2array
# from histo_utilities import create_TH2D, create_TH1D
from cebefo_style import cebefo_style
import os, sys, re, glob
import matplotlib.pyplot as plt

rt.gErrorIgnoreLevel = 4000

trigger = 'HT'
# trigger = 'TOF'

if len(sys.argv)>1:
    trigger = sys.argv[1]
if not (trigger in ['HT', 'TOF']):
    print 'Available triggers: HT, TOF'
    raise


cebefo_style()
c_light = 2.99792458E8 #m/s
lumi = 12.3 # fb-1
donotdelete = []

colors = [1, 2, 4, 6 ,8, 5, 13, 46]

if trigger == 'HT':
    trigger_selection = '(max_vtx_SumPT > 350)'
elif trigger == 'TOF':
    trigger_selection = '(mh>10 && same_vtx_idx && max_vtx_SumPT > 150)'

save_dir = '_fig/' + trigger + 'trigger'
print 'Running with trigger -->', trigger

min_trk_pt = [100., 50.]
max_dm = 0.25

mass_lower_bound = 50
mass_upper_bound = 600
binning_pt = [50, 30, 800]
arr = np.concatenate((np.array([-1, 0]), np.logspace(-2, 1, 50)))
binning_dm = [arr.shape[0]-1, arr]

binning_mass = [(mass_upper_bound-mass_lower_bound)/10, mass_lower_bound, mass_upper_bound]

labels = []
production = {} # tracks in the mass spectrum per fb-1 of integrated lumi
bins_pt_hat = []

file_path = glob.glob('_root/flat_evt_copies/pp2HardQCD_PU140_pTHat*')
Ttrees = {}
Tfiles = {}
for path in file_path:
    out = re.search(r'at[0-9]+-[0-9,I, N, F]+', os.path.basename(path[:-17]))
    label = out.group(0)[2:]
    labels.append(label)
    print label

    bins_pt_hat.append(float(label[:3]))
    if label[-3:] == 'INF':
        bins_pt_hat.append(2*bins_pt_hat[-1])


    Tfiles[label] = rt.TFile(path, 'READ')
    Ttrees[label] = Tfiles[label].Get('T')

    production[label] = [Tfiles[label].Get('XSec')[0], Tfiles[label].Get('XSec')[1], Tfiles[label].Get('NEvts')[0]]

lab_int = [int(l[:3])for l in labels]
labels = np.array(labels)
labels = list(labels[np.argsort(lab_int)])

h_mass_cat1 = {}
h_mass_cat2 = {}
h_dm_pth = {}

leg = rt.TLegend(0.7, 0.65, 0.95, 0.98)
weight_sum = 0

for i, l in enumerate(labels):
    t = Ttrees[l]

    evt_weight = lumi * production[l][0] / production[l][2]
    weight_sum += evt_weight


    h_dm_pth[l] = rt.TH2D('h_dm_pth_'+l, l, binning_pt[0], binning_pt[1], binning_pt[2], binning_dm[0], binning_dm[1])
    t.Project('h_dm_pth_'+l, '(mh/ml - 1)*(ml>0) - (ml<0) : pth', trigger_selection)
    h_dm_pth[l].Scale(evt_weight)
    h_dm_pth[l].evt_weight = evt_weight


    dm_cut = '(mh/ml - 1 > {} || ml < 0)'.format(max_dm)
    sel_cat1 = trigger_selection + ' && ' + dm_cut
    sel_cat1 += ' && (pth > {})'.format(min_trk_pt[0])

    h_mass_cat1[l] = rt.TH1D('h_mass_cat1_'+l, l, binning_mass[0], binning_mass[1], binning_mass[2])
    t.Project('h_mass_cat1_'+l, 'mh', sel_cat1)
    h_mass_cat1[l].Scale(evt_weight)


    dm_cut = '(mh/ml - 1 < {} && ml > 0)'.format(max_dm)
    sel_cat2 = trigger_selection + ' && ' + dm_cut
    sel_cat2 += ' && (pth > {})'.format(min_trk_pt[1])

    h_mass_cat2[l] = rt.TH1D('h_mass_cat2_'+l, l, binning_mass[0], binning_mass[1], binning_mass[2])
    t.Project('h_mass_cat2_'+l, '(mh+ml)/2', sel_cat2)
    h_mass_cat2[l].Scale(evt_weight)

line = rt.TLine()
line.SetLineColor(1)
line.SetLineWidth(3)
line.SetLineStyle(9)

c_dm_pth = rt.TCanvas("c_dm_pth", "c_dm_pth", 800, 600)
c_dm_pth.dnd = []
h_dm_pth_tot = rt.TH2D('h_dm_pth', 'dm VS pth', binning_pt[0], binning_pt[1], binning_pt[2], binning_dm[0], binning_dm[1])
h_dm_pth_tot.SetXTitle('p_{T}^{h} [GeV]')
h_dm_pth_tot.SetYTitle('m_{h}/m_{l} - 1')
for i, l in enumerate(labels):
    h_dm_pth_tot.Add(h_dm_pth[l])
    # h_aux = h_dm_pth_tot.Clone('hauxpt_'+l)
    # h_aux.SetLineWidth(2)
    # h_aux.SetLineColor(colors[-i-1])
    # h_aux.SetMarkerColor(colors[-i-1])
    # h_aux.SetStats(0)
    # c_dm_pth.cd()
    # h_aux.Draw("same")
    # c_dm_pth.dnd.append(h_aux)
# c_dm_pth.dnd[0].GetYaxis().SetRangeUser(1., h_dm_pth_tot.GetMaximum())
# h_dm_pth_tot.SetStats(0)
h_dm_pth_tot.Draw('colz')
# leg.Draw()
# line.DrawLine(min_trk_pt[0], 1., min_trk_pt[0], h_dm_pth_tot.GetMaximum())
rt.gPad.SetRightMargin(0.99)
c_dm_pth.SetLogz()
c_dm_pth.SaveAs(save_dir + '/dm_pth_cat1.png')
c_dm_pth.Update()


c_mass_cat1 = rt.TCanvas("c_mass_cat1", "c_mass_cat1", 800, 600)
c_mass_cat1.dnd = []
h_mass_cat1_tot = rt.TH1D('h_mass_cat1', 'mass cat1', binning_mass[0], binning_mass[1], binning_mass[2])
h_mass_cat1_tot.SetXTitle('m_{h} [GeV]')
h_mass_cat1_tot.SetYTitle('Events / {:.1f} GeV'.format(h_mass_cat1_tot.GetBinWidth(1)))
for i, l in enumerate(labels):
    h_mass_cat1_tot.Add(h_mass_cat1[l])
    h_aux = h_mass_cat1_tot.Clone('hauxpt_'+l)
    h_aux.SetLineWidth(2)
    h_aux.SetLineColor(colors[-i-1])
    h_aux.SetMarkerColor(colors[-i-1])
    h_aux.SetStats(0)
    leg.AddEntry(h_aux, '#hat{p}_{T} < ' + l[-3:] + ' GeV', 'lep')
    c_mass_cat1.cd()
    h_aux.Draw("same")
    c_mass_cat1.dnd.append(h_aux)
c_mass_cat1.dnd[0].GetYaxis().SetRangeUser(1, h_mass_cat1_tot.GetMaximum())

# f_fit_sum = rt.TF1('f_fit_sum', '[0] * exp(-[1]*x)', mass_lower_bound, binning_mass[2])
# f_fit_sum = rt.TF1('f_fit_sum', '[0] * pow(x, -[1])', mass_lower_bound, binning_mass[2])
f_fit_sum = rt.TF1('f_fit_sum', '[0] * exp (-[1]*x)/ x ', mass_lower_bound, binning_mass[2])
f_fit_sum.SetParameter(0, 2*h_mass_cat1_tot.GetBinContent(0))
f_fit_sum.SetParameter(1, 1./65.)
res = h_mass_cat1_tot.Fit('f_fit_sum','SIL0R+')
chi2 = 0
for i in range(h_mass_cat1_tot.GetNbinsX()):
    j = i+1
    if h_mass_cat1_tot.GetBinError(j) != 0:
        chi2 += (h_mass_cat1_tot.GetBinContent(j) - f_fit_sum.Eval(h_mass_cat1_tot.GetBinCenter(j)))**2 / h_mass_cat1_tot.GetBinError(j)**2
    # print j, chi2
f_fit_sum.SetLineStyle(9)
f_fit_sum.SetLineColor(1)
f_fit_sum.Draw('SameL')
note = rt.TLatex()
note.SetTextSize(0.04)
# note.DrawLatexNDC(0.3, 0.88, '#chi^{{2}}/ndf = {:.1f}/{}'.format(chi2, res.Get().Ndf()))
Mstar1 = 1./f_fit_sum.GetParameter(1)
DMstar1 = f_fit_sum.GetParError(1)/Mstar1**2
note.DrawLatexNDC(0.3, 0.90, 'M* = {:.1f} GeV'.format(Mstar1))
N_bkg_evts1 = h_mass_cat1_tot.Integral()
note.DrawLatexNDC(0.3, 0.85, 'Events (L = {:.1f} fb^{{-1}}) = {:1.2e}'.format(lumi, h_mass_cat1_tot.Integral()))

leg.Draw()
rt.gPad.SetRightMargin(0.99)
c_mass_cat1.SetLogy()
c_mass_cat1.SaveAs(save_dir + '/mass_cat1.png')
c_mass_cat1.SaveAs(save_dir + '/mass_cat1.root')
c_mass_cat1.Update()



c_mass_cat2 = rt.TCanvas("c_mass_cat2", "c_mass_cat2", 800, 600)
c_mass_cat2.dnd = []
h_mass_cat2_tot = rt.TH1D('h_mass_cat2', 'mass cat2', binning_mass[0], binning_mass[1], binning_mass[2])
h_mass_cat2_tot.SetXTitle('(m_{h} + m_{l})/2 [GeV]')
h_mass_cat2_tot.SetYTitle('Events / {:.1f} GeV'.format(h_mass_cat2_tot.GetBinWidth(1)))
for i, l in enumerate(labels):
    h_mass_cat2_tot.Add(h_mass_cat2[l])
    h_aux = h_mass_cat2_tot.Clone('hauxpt_'+l)
    h_aux.SetLineWidth(2)
    h_aux.SetLineColor(colors[-i-1])
    h_aux.SetMarkerColor(colors[-i-1])
    h_aux.SetStats(0)
    c_mass_cat2.cd()
    h_aux.Draw("same")
    c_mass_cat2.dnd.append(h_aux)
c_mass_cat2.dnd[0].GetYaxis().SetRangeUser(1, h_mass_cat2_tot.GetMaximum())

f_fit_sum = rt.TF1('f_fit_sum', '[0] * exp(-[1]*x)', mass_lower_bound, binning_mass[2])
# f_fit_sum = rt.TF1('f_fit_sum', '[0] * pow(x, -[1])', mass_lower_bound, binning_mass[2])
# f_fit_sum = rt.TF1('f_fit_sum', '[0] * exp (-[1]*x)/ x ', mass_lower_bound, binning_mass[2])
f_fit_sum.SetParameter(0, h_mass_cat2_tot.Integral()/35.)
f_fit_sum.SetParameter(1, 1./35.)
res = h_mass_cat2_tot.Fit('f_fit_sum','SIL0R+')
chi2 = 0
for i in range(h_mass_cat2_tot.GetNbinsX()):
    j = i+1
    if h_mass_cat2_tot.GetBinError(j) != 0:
        chi2 += (h_mass_cat2_tot.GetBinContent(j) - f_fit_sum.Eval(h_mass_cat2_tot.GetBinCenter(j)))**2 / h_mass_cat2_tot.GetBinError(j)**2
    # print j, chi2
f_fit_sum.SetLineStyle(9)
f_fit_sum.SetLineColor(1)
f_fit_sum.Draw('SameL')
note = rt.TLatex()
note.SetTextSize(0.04)
note.DrawLatexNDC(0.3, 0.88, '#chi^{{2}}/ndf = {:.1f}/{}'.format(chi2, res.Get().Ndf()))
# note.DrawLatexNDC(0.3, 0.88, '#chi^{{2}}/ndf = {:.1f}/{} ({:.2f})'.format(res.Get().Chi2(), res.Get().Ndf(), res.Get().Prob()))
Mstar2 = 1./f_fit_sum.GetParameter(1)
DMstar2 = f_fit_sum.GetParError(1)/Mstar2**2
note.DrawLatexNDC(0.3, 0.83, 'M* = {:.1f} GeV'.format(Mstar2))
N_bkg_evts2 = h_mass_cat2_tot.Integral()
note.DrawLatexNDC(0.3, 0.78, 'Events (L = {:.1f} fb^{{-1}}) = {:1.2e}'.format(lumi, h_mass_cat2_tot.Integral()))

leg.Draw()
rt.gPad.SetRightMargin(0.99)
c_mass_cat2.SetLogy()
c_mass_cat2.SaveAs(save_dir + '/mass_cat2.png')
c_mass_cat2.Update()


f_pars = open(save_dir + '/bkg_parametrization.txt', 'w')
f_pars.write('double Lumi_Nbkg = {:1.2e}; //[pb-1]\n'.format(lumi*1e3))
f_pars.write('double Nevs_bkg_{}[] = {{{:1.2e}, {:1.2e}}};\n'.format(trigger, N_bkg_evts1, N_bkg_evts2))
f_pars.write('double Mstar_bkg_{}[] = {{{:.1f}, {:.1f}}};\n'.format(trigger, Mstar1, Mstar2))
f_pars.write('double DMstar_bkg_{}[] = {{{:.1g}, {:.1g}}};\n'.format(trigger, DMstar1, DMstar2))

import numpy as np
import ROOT as rt
from root_numpy import root2array
from histo_utilities import create_TH2D, create_TH1D
from cebefo_style import cebefo_style
import os, sys, re, glob
import matplotlib.pyplot as plt

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
    out = re.search(r'_[0-9]+k', path)
    out = float(out.group(0)[1:-1])*1000
    print 'Events gen = {:d}k'.format(int(out/1000))

    return (xsec, err, out)


def compute_Eff(p, t, alpha = 0.05):
    e = p/t
    if e == 0:
        de = 1 - np.power(alpha, 1/t)
    else:
        de = np.sqrt(p*(1-e))/t
    return e,de


cebefo_style()
c_light = 2.99792458E8 #m/s
Luminosity = 1000 # fb-1
donotdelete = []

min_trk_pt = 100
max_dm = 0.25
min_vtx_pt = 350
mass_lower_bound = 50
binning_vtx = [100, 0, 3000, 100, 0, 800]

binning_mass = [(600-mass_lower_bound)/10, mass_lower_bound, 600.]


file_path = list(sys.argv[1:])
labels = []
production = {} # tracks in the mass spectrum per fb-1 of integrated lumi
bins_pt_hat = []

h_mass = {}

c_summary = rt.TCanvas('c_summarybins', 'c_summarybins', 1600, 1600)
c_summary.Divide(2,len(file_path))

weight_sum = 0

for i, path in enumerate(file_path):
    name = os.path.basename(path[:-14])
    print name

    out = re.search(r'at[0-9]+-[0-9,I, N, F]+', name)
    label = out.group(0)[2:]
    labels.append(label)

    bins_pt_hat.append(float(label[:3]))
    if label[-3:] == 'INF':
        bins_pt_hat.append(2*bins_pt_hat[-1])


    t = root2array(path)
    dm = t['mh']/t['ml'] - 1
    mass_sel = np.logical_or(dm > max_dm, t['ml'] < 0)

    arr_pt_vtxpt = np.column_stack((t['vtx_SumPT'][mass_sel], t['pth'][mass_sel]))
    h_vtx = create_TH2D(arr_pt_vtxpt, 'h_pt_vtxpt_'+name, name, binning=binning_vtx, axis_title=['#sum p_{T} [GeV]','p_{T}^{H} [GeV]'])


    kin_sel = np.logical_and(np.abs(t['pth']) >= min_trk_pt, np.abs(t['vtx_SumPT']) >= min_vtx_pt)
    kin_sel = np.logical_and(kin_sel, mass_sel)
    kin_sel = np.logical_and(kin_sel, t['mh'] > mass_lower_bound)

    arr_mh = t['mh'][kin_sel]
    h_mass[label] = create_TH1D(arr_mh, 'h_mass_'+name, name, binning=binning_mass, axis_title=['(m_{H} [GeV]','Events / {:.1} GeV'.format(float(binning_mass[2]-binning_mass[1])/binning_mass[0])])


    '''  ------------------ Drawing --------------- '''
    line = rt.TLine()
    line.SetLineColor(2)
    line.SetLineWidth(2)
    line.SetLineStyle(7)

    c_summary.cd(2*i + 1)
    h_vtx.DrawCopy('colz')
    line.DrawLine(min_vtx_pt, min_trk_pt, min_vtx_pt, binning_vtx[5])
    line.DrawLine(min_vtx_pt, min_trk_pt, binning_vtx[2], min_trk_pt)

    c_summary.cd(2*i + 2)
    h_mass[label].Draw('E1')
    rt.gPad.SetLogy()

    # Getting the xsec
    production[label] = GetXsec(path)
    print 'Passed: {:d}'.format(int(h_mass[label].GetEntries()))

    #Computing the weights
    h_mass[label].weight = h_mass[label].GetEntries() * production[label][0] / production[label][2]
    weight_sum += h_mass[label].weight
    print 'Weight: {:.1f}'.format(h_mass[label].weight)
    print '\n'

c_summary.Update()
c_summary.SaveAs('/Users/olmo/Desktop/bkg_summary.png')

# Draw a weighted sum of the bkg
c_sum = rt.TCanvas("c_sum", "c_sum", 800, 600)
c_sum.donotdelete = []
h_stack = create_TH1D([], 'h_mass_sum', 'Mass spectrum of HardQCD weighted', binning=binning_mass, axis_title=['Mass reco [GeV]', 'Norm. entries / {:.1f} GeV'.format(binning_mass[2]/binning_mass[0])])
h_stack.Sumw2()
leg = rt.TLegend(0.7, 0.65, 0.95, 0.98)
colors = [1, 2, 4, 6 ,8, 5]
for i, l in enumerate(labels):
    if h_mass[l].GetEntries() == 0:
        continue
    weight = h_mass[l].weight / (h_mass[l].GetEntries() * weight_sum)
    h_stack.Add(h_mass[l], weight)
    h_tmp = h_stack.Clone("drawn_tmp_"+l)
    h_tmp.SetTitle('')
    h_tmp.GetYaxis().SetRangeUser(1.e-6, 1)
    h_tmp.GetYaxis().SetTitleOffset(1.15)
    h_tmp.SetStats(0)
    h_tmp.SetLineColor(colors[i])
    h_tmp.SetMarkerColor(colors[i])
    leg.AddEntry(h_tmp, '#hat{p}_{T} < ' + l[-3:] + ' GeV', 'lep')
    if i == 0:
        h_tmp.Draw("E1")
    else:
        h_tmp.Draw("SAME E1")

    c_sum.donotdelete.append(h_tmp)

# f_fit_sum = rt.TF1('f_fit_sum', '[0] * exp(-[1]*x)', mass_lower_bound, binning_mass[2])
# f_fit_sum = rt.TF1('f_fit_sum', '[0] * pow(x, -[1])', mass_lower_bound, binning_mass[2])
f_fit_sum = rt.TF1('f_fit_sum', '[0] * exp (-[1]*x)/ x ', mass_lower_bound, binning_mass[2])
f_fit_sum.SetParameter(0, mass_lower_bound)
f_fit_sum.SetParameter(1, 1./mass_lower_bound)
res = h_stack.Fit('f_fit_sum','SL0R+')
f_fit_sum.Draw('SameL')
note = rt.TLatex()
note.DrawLatexNDC(0.2, 0.88, '#chi^{{2}}/ndf = {:.1f}/{} ({:.2f})'.format(res.Get().Chi2(), res.Get().Ndf(), res.Get().Prob()))
leg.Draw()
c_sum.SetLogy()
rt.gPad.SetLeftMargin(.13)
rt.gPad.SetRightMargin(.05)
rt.gPad.SetBottomMargin(.14)
rt.gPad.SetTopMargin(.02)
c_sum.Update()
print 'Normalized bkg fit parameters', f_fit_sum.GetParameter(0)/f_fit_sum.Integral(mass_lower_bound, 1e5), f_fit_sum.GetParameter(1)
print 'M* = ', 1./f_fit_sum.GetParameter(1)
c_sum.SaveAs('/Users/olmo/Desktop/bkg_sum.png')


bins_pt_hat = np.sort(np.array(bins_pt_hat))
h_tot = rt.TH1F('h_tot', 'h_tot', bins_pt_hat.shape[0]-1, bins_pt_hat)
h_pass = rt.TH1F('h_pass', 'h_pass', bins_pt_hat.shape[0]-1, bins_pt_hat)

xsec = np.zeros(bins_pt_hat.shape[0]-1)

for l in labels:
    bin_number = np.argmax(bins_pt_hat == float(l[:3])) + 1

    h_tot.SetBinContent(bin_number, production[l][2])
    h_pass.SetBinContent(bin_number, h_mass[l].GetEntries())

    xsec[bin_number-1] = production[l][0]

eff_calc = rt.TEfficiency(h_pass, h_tot)

gr = eff_calc.CreateGraph()
tot_evt = [0, 0]
print 'p_T(Hat)    N_exp'
for i in range(gr.GetN()):
    scale = xsec[i] * Luminosity

    x = rt.Double(0)
    y = rt.Double(0)

    gr.GetPoint(i, x, y)
    gr.SetPoint(i, x, y*scale)
    print '{:.1f}    {:1.2e}'.format(x, y*scale)
    tot_evt[0] +=  y*scale
    tot_evt[1] +=  y*xsec[i]*300
    gr.SetPointEYhigh(i, gr.GetErrorYhigh(i) * scale)
    gr.SetPointEYlow(i, gr.GetErrorYlow(i) * scale)

print "Tot bkg evt @{} fb-1 = {:.2e}".format(Luminosity, tot_evt[0])
print "Tot bkg evt @300 fb-1 = {:.2e}".format(tot_evt[1])
print "Tot bkg evt @12 fb-1 = {:.2e}".format(12.*tot_evt[1]/300.)

c_prod = rt.TCanvas('c_prod', 'c_prod', 800, 600)
gr.SetFillColor(8)
gr.SetFillStyle(3001)
gr.GetXaxis().SetTitle('#hat{p}_{T} [GeV]')
gr.GetYaxis().SetTitle('Number of events')
gr.GetYaxis().SetLimits(1., 1.e9)
gr.GetXaxis().SetRangeUser(np.min(bins_pt_hat), np.max(bins_pt_hat))
gr.Draw('A20')
gr.Draw('P')

gr_MC = rt.TGraphErrors(gr.GetN())
for l in labels:
    bin_number = np.argmax(bins_pt_hat == float(l[:3])) + 1
    m_min = float(l[:3])
    if l[-3:] == 'INF':
        m_max = 2*m_min
    else:
        m_max = float(l[-3:])

    x = 0.5*(m_max + m_min)
    gr_MC.SetPoint(bin_number, x, h_mass[l].GetEntries())
    gr_MC.SetPointError(bin_number, 0.5*(m_max - m_min), 0)


gr_MC.SetLineColor(2);
gr_MC.SetLineWidth(3);
gr_MC.SetMarkerColor(2);
gr_MC.Draw('P')

c_prod.SetLogy()
c_prod.SetLogx()
c_prod.Update()
c_prod.SaveAs('/Users/olmo/Desktop/production_stats.root')

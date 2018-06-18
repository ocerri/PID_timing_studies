import numpy as np
import ROOT as rt
from array import array
# from root_numpy import root2array
# from histo_utilities import create_TH2D, create_TH1D
from cebefo_style import cebefo_style
import os, sys, re, glob
import matplotlib.pyplot as plt

rt.gErrorIgnoreLevel = 4000

trigger = 'HT'
# trigger = 'TOF'


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

binning_pt = [50, 30, 800]
arr = np.concatenate((np.array([-0.5, -0.3]), np.linspace(0.0, 0.6, 20)))
binning_dm = [arr.shape[0]-1, arr]

stopmass = []

file_path = glob.glob('_root/flat_evt_copies/pp2RHad_PU140_M*')
Ttrees = {}
Tfiles = {}
N_evts = {}
for path in file_path:
    name = os.path.basename(path)
    out = re.search('_M[0-9]+_', name)
    out = out.group(0)
    m = int(out[2:-1])
    if m <= 50 or m >= 3000:
        continue
    stopmass.append(m)

    Tfiles[m] = rt.TFile(path, 'READ')
    Ttrees[m] = Tfiles[m].Get('T')

    out = re.search(r'_[0-9]+k', name)
    N_evts[m] = float(out.group(0)[1:-1])*1000

stopmass =  list(np.sort(stopmass))

pars = {1: [], 2: []}

for i, m in enumerate(stopmass):
    t = Ttrees[m]
    l = str(m)

    c_out =  rt.TCanvas('c_out'+l, 'c_out'+l, 1200, 600)
    c_out.Divide(3, 1)

    h_dm_pth = rt.TH2D('h_dm_pth_'+l, l, binning_pt[0], binning_pt[1], 3*m, binning_dm[0], binning_dm[1])
    t.Project('h_dm_pth_'+l, '(mh/ml - 1)*(ml>0) - 0.5*(ml<0) : pth', trigger_selection)
    h_dm_pth.Scale(1./h_dm_pth.Integral())
    c_out.cd(1)
    h_dm_pth.SetXTitle('p_{T}^{h} [GeV]')
    h_dm_pth.SetYTitle('m_{h}/m_{l} - 1')
    h_dm_pth.Draw('colz')
    rt.gPad.SetLogz()


    # -------------------------- Gauss Double exp TFormula ---------------------
    '''
    [0] = Norm
    [1] = aL
    [2] = aR
    [3] = mu
    [4] = sigma
    '''
    # Norm = 'sqrt(3.14159/2.)*[4]*( erf([2]/sqrt(2)) - erf(-[1]/sqrt(2)) ) '
    # Norm += '+ exp(-0.5*[1]*[1])*[4]/[1] + exp(-0.5*[2]*[2])*[4]/[2]'
    Norm = '1'

    t_str = '((x-[3])/[4])'
    left_fun = '(TMath::Sign(0.5,-{}-[1])+0.5)*exp(0.5*[1]*[1])*exp([1]*{})'.format(t_str,t_str)
    right_fun = '(TMath::Sign(0.5,{}-[2])+0.5)*exp(0.5*[2]*[2])*exp(-[2]*{})'.format(t_str, t_str)
    central_fun = '(TMath::Sign(0.5,{}+[1])+0.5)*(TMath::Sign(0.5,[2]-{})+0.5)*exp(-0.5*{}*{})'.format(t_str, t_str, t_str, t_str)
    TF_gausDE = '[0]*(' + Norm + ')*(' + left_fun + " + " + central_fun + ' + ' + right_fun + ')'
    # --------------------------------------------------------------------------

    # ----------------- Category #1 --------------------------------------------
    binning_mass = [80, m/1.8, m*1.6]

    dm_cut = '(mh/ml - 1 > {} || ml < 0)'.format(max_dm)
    sel_cat1 = trigger_selection + ' && ' + dm_cut
    sel_cat1 += ' && (pth > {})'.format(min_trk_pt[0])

    h_mass_cat1 = rt.TH1D('h_mass_cat1_'+l, 'Category 1, ' + l, binning_mass[0], binning_mass[1], binning_mass[2])
    t.Project('h_mass_cat1_'+l, 'mh', sel_cat1)
    eff1 = h_mass_cat1.GetEntries() / float(N_evts[m])
    eff1_err = np.sqrt(eff1 * (1 - eff1) / float(N_evts[m]) )

    h_mass_cat1.Scale(1./h_mass_cat1.Integral())
    q = np.zeros(2)
    h_mass_cat1.GetQuantiles(2, q, np.array([0.02, 0.98]))
    f_fit1 = rt.TF1('f_fit1_'+l, TF_gausDE, q[0] , q[1])
    f_fit1.SetParNames('K', '#alpha_{L}', '#alpha_{R}' , '#mu', '#sigma')

    integral = h_mass_cat1.Integral()*h_mass_cat1.GetBinWidth(1)
    # print integral
    # integral /= np.sqrt(2*np.pi)*h_mass_cat1.GetRMS()
    f_fit1.SetParameter(0, 0.15) #Norm
    f_fit1.SetParLimits(0,1e-3, 50)
    # f_fit1.FixParameter(0, 0.15) #Norm

    f_fit1.SetParameter(1, 1.5) #left alpha
    f_fit1.SetParLimits(1, 0.5, 5)
    # f_fit1.FixParameter(1, 1.5) #left alpha

    f_fit1.SetParameter(2, 1.5) #right alpha
    f_fit1.SetParLimits(2, 0.5, 5)
    # f_fit1.FixParameter(2, 1.5) #right alpha

    f_fit1.SetParameter(3,h_mass_cat1.GetMean()) #mean
    f_fit1.SetParLimits(3,binning_mass[1], binning_mass[2])
    # f_fit1.FixParameter(3,h_mass_cat1.GetMean()) #mean

    f_fit1.SetParameter(4, h_mass_cat1.GetRMS()*0.7) #sigma
    f_fit1.SetParLimits(4, 1, 1000)
    # f_fit1.FixParameter(4, h_mass_cat1.GetRMS()*0.7) #right alpha

    # print 'Integral:', f_fit1.Integral(binning_mass[0], binning_mass[1])

    fit_results = h_mass_cat1.Fit(f_fit1, 'IL0SQR+')

    pars[1].append([m,
                 float(f_fit1.GetParameter(1)),
                 float(f_fit1.GetParError(1)),
                 float(f_fit1.GetParameter(2)),
                 float(f_fit1.GetParError(2)),
                 float(f_fit1.GetParameter(3)),
                 float(f_fit1.GetParError(3)),
                 float(f_fit1.GetParameter(4)),
                 float(f_fit1.GetParError(4)),
                 fit_results.Chi2(),
                 int(fit_results.Ndf()),
                 eff1,
                 eff1_err
                 ])

    donotdelete.append(f_fit1)
    f_fit1.SetLineColor(2)
    f_fit1.SetLineWidth(2)
    f_fit1.SetLineStyle(9)
    c_out.cd(2)
    h_mass_cat1.SetXTitle('m_{h} [GeV]')
    h_mass_cat1.SetYTitle('Prob. / {:.1f} GeV'.format(h_mass_cat1.GetBinWidth(1)))
    h_mass_cat1.Draw()
    f_fit1.Draw("SAME")
    # rt.gPad.SetLogy()

    # ----------------- Category #2 --------------------------------------------
    binning_mass = [50, m/1.4, m*1.3]

    dm_cut = '(mh/ml - 1 < {} && ml > 0)'.format(max_dm)
    sel_cat2 = trigger_selection + ' && ' + dm_cut
    sel_cat2 += ' && (pth > {})'.format(min_trk_pt[1])

    h_mass_cat2 = rt.TH1D('h_mass_cat2_'+l, 'Category 2, ' + l, binning_mass[0], binning_mass[1], binning_mass[2])
    t.Project('h_mass_cat2_'+l, '(mh+ml)/2', sel_cat2)
    eff2 = h_mass_cat2.GetEntries() / float(N_evts[m])
    eff2_err = np.sqrt(eff2 * (1 - eff2) / float(N_evts[m]) )
    h_mass_cat2.Scale(1./h_mass_cat2.Integral())

    q = np.zeros(2)
    h_mass_cat2.GetQuantiles(2, q, np.array([0.03, 0.98]))
    f_fit = rt.TF1('f_fit1_'+l, TF_gausDE, q[0] , q[1])
    f_fit.SetParNames('K', '#alpha_{L}', '#alpha_{R}' , '#mu', '#sigma')

    integral = h_mass_cat2.Integral()*h_mass_cat2.GetBinWidth(1)
    integral /= np.sqrt(2*np.pi)*h_mass_cat2.GetRMS()
    f_fit.SetParameter(0, integral) #Norm
    f_fit.SetParLimits(0, 1e-4, 100)

    f_fit.SetParameter(1,1.7) #left alpha
    f_fit.SetParLimits(1, 0.5, 5)
    # f_fit.FixParameter(1,5) #left alpha

    f_fit.SetParameter(2, 1.7) #right alpha
    f_fit.SetParLimits(2, 0.5, 5)
    # f_fit.FixParameter(2,5) #right alpha

    f_fit.SetParameter(3,h_mass_cat2.GetMean()) #mean
    f_fit.SetParLimits(3,binning_mass[1], binning_mass[2])

    f_fit.SetParameter(4, h_mass_cat2.GetRMS()*0.8) #sigma
    f_fit.SetParLimits(4, 1, 1000)

    fit_results = h_mass_cat2.Fit(f_fit, 'IL0SQR+')

    pars[2].append([m,
                 float(f_fit.GetParameter(1)),
                 float(f_fit.GetParError(1)),
                 float(f_fit.GetParameter(2)),
                 float(f_fit.GetParError(2)),
                 float(f_fit.GetParameter(3)),
                 float(f_fit.GetParError(3)),
                 float(f_fit.GetParameter(4)),
                 float(f_fit.GetParError(4)),
                 fit_results.Chi2(),
                 fit_results.Ndf(),
                 eff2,
                 eff2_err
                 ])

    donotdelete.append(f_fit)
    f_fit.SetLineColor(2)
    f_fit.SetLineWidth(2)
    f_fit.SetLineStyle(9)
    c_out.cd(3)
    h_mass_cat2.SetXTitle('(m_{h} + m_{l})/2 [GeV]')
    h_mass_cat2.SetYTitle('Prob. / {:.1f} GeV'.format(h_mass_cat2.GetBinWidth(1)))
    h_mass_cat2.Draw()
    f_fit.Draw("SAME")
    # rt.gPad.SetLogy()

    # ------------------ Save stuff --------------------------
    c_out.Update()
    donotdelete.append([c_out, h_dm_pth, h_mass_cat1, h_mass_cat2])
    c_out.SaveAs(save_dir + '/RHad_mass{}_spectrum.png'.format(m))

    if m == 500:
        print m
        c_out2 = rt.TCanvas('c_out2_'+l, 'c_out2_'+l, 1200, 600)
        c_out2.Divide(2, 1)
        c_out2.cd(1)
        h_mass_cat1.Draw()
        f_fit1.Draw("SAME")
        c_out2.cd(2)
        h_mass_cat2.Draw()
        f_fit.Draw("SAME")
        c_out2.SaveAs(save_dir + '/RHad_mass{}_spectrum.root'.format(m))



pars[1] = np.array(pars[1])
pars[2] = np.array(pars[2])
pars_names = ['M_st', 'aL', 'aL_err', 'aR', 'aR_err', 'mu', 'mu_err', 'sigma', 'sigma_err', 'chi2', 'NDF', 'eff']

# -------------------- Dump the parameters in a file --------------
f_pars = open(save_dir + '/signal_parametrization.txt', 'w')
for i in [1,2]:
    f_pars.write('Category {}\n'.format(i))
    for j, n in enumerate(pars_names):
        if 'err' in n:
            continue
        s = ''
        for x in pars[i][:, j]:
            s += '{:.3f}, '.format(x)
        f_pars.write('{}_{}{}'.format(n, trigger, i) + ' = {' + s[:-2] + '};\n')
    f_pars.write('\n\n'.format(i))
f_pars.close()


fig = plt.figure()
plt.errorbar(pars[1][:,0], pars[1][:,10], yerr=np.sqrt(pars[1][:,10]), fmt='r.', label = 'NDF Cat 1')
plt.plot(pars[1][:,0], pars[1][:,9], 'rx')
plt.errorbar(pars[2][:,0], pars[2][:,10], yerr=np.sqrt(pars[2][:,10]), fmt='b.', label = 'NDF Cat 2')
plt.plot(pars[2][:,0], pars[2][:,9], 'bx')
plt.ylabel(r'$\chi^{2}$')
plt.xlabel(r'$M_{stop}$ [GeV]')
# plt.yscale("log")
plt.grid()
plt.savefig(save_dir + '/Chi2.png')
plt.legend(loc='best')
fig.show()

c = rt.TCanvas('c', 'c', 800, 600)
c.Divide(1, 2)

Npts = pars[1].shape[0]

c.cd(1)
g_s1 = rt.TGraph(Npts)
# g_s1 = rt.TGraphErrors(Npts)
for i in range(Npts):
    g_s1.SetPoint(i, pars[1][i, 0], pars[1][i, 7])
    # g_s1.SetPointError(i, 0, pars[1][i, 8])

g_s2 = rt.TGraph(Npts)
# g_s2 = rt.TGraphErrors(Npts)
for i in range(Npts):
    g_s2.SetPoint(i, pars[2][i, 0], pars[2][i, 7])
    # g_s2.SetPointError(i, 0, pars[2][i, 8])

g_s1.GetXaxis().SetTitle('#tilde{t}_1 mass [GeV]')
g_s1.GetYaxis().SetTitle('GaussDoubleExp #sigma [GeV]')
g_s1.GetYaxis().SetRangeUser(1, 300)
g_s1.SetTitle('')
g_s1.SetLineColor(2)
g_s1.SetMarkerColor(2)
g_s1.SetMarkerStyle(4)
g_s1.Draw('APE1C')
g_s2.SetLineColor(4)
g_s2.SetMarkerColor(4)
g_s2.SetMarkerStyle(4)
rt.gPad.SetLogy()
g_s2.Draw('PE1C')


c.cd(2)
g_e1 = rt.TGraphErrors(Npts)
for i in range(Npts):
    g_e1.SetPoint(i, pars[1][i, 0], pars[1][i, 11])
    g_e1.SetPointError(i, 0, pars[1][i, 12])

g_e2 = rt.TGraphErrors(Npts)
for i in range(Npts):
    g_e2.SetPoint(i, pars[2][i, 0], pars[2][i, 11])
    g_e2.SetPointError(i, 0, pars[2][i, 12])

g_e1.GetXaxis().SetTitle('#tilde{t}_1 mass [GeV]')
g_e1.GetYaxis().SetTitle('Efficiency')
g_e1.GetYaxis().SetRangeUser(0.005, 1)
g_e1.SetTitle('')
g_e1.SetLineColor(2)
g_e1.SetMarkerColor(2)
g_e1.SetMarkerStyle(4)
g_e1.Draw('APE1C')
g_e2.SetLineColor(4)
g_e2.SetMarkerColor(4)
g_e2.SetMarkerStyle(4)
g_e2.Draw('PE1C')
rt.gPad.SetLogy()

leg = rt.TLegend(0.8, 0.2, 0.98, 0.48)
leg.AddEntry(g_e1, 'Cat 1', 'lep')
leg.AddEntry(g_e2, 'Cat 2', 'lep')
leg.Draw()

c.Update()
c.SaveAs(save_dir + '/SigmaEff.root')

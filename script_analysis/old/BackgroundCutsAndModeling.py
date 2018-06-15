import numpy as np
import ROOT as rt
from root_numpy import root2array
from histo_utilities import create_TH2D, create_TH1D
from cebefo_style import cebefo_style
import os, sys, re, glob
import matplotlib.pyplot as plt
# plt.style.use('olmo')

cebefo_style()
c_light = 2.99792458E8 #m/s
donotdelete = []

min_trk_pt = 100
min_vtx_pt = 370

file_path = list(sys.argv[1:])
print len(sys.argv)-1

pars = []
xsec = [2.555e-2, 7.113e-05]
# log_path =  path[:-4] + 'log'
# cmd = 'more ' + log_path + ' | grep "| sum"'
# out = os.system(cmd)
# print 'Ciao: ', out

tree = {}
h_pt = {}
binning = [100, 0, 3000, 100, 0, 1500]
axis_title = ['reco vtx sum p_{T} [GeV]', 'reco p_{T} [GeV]']

for path in file_path:
    name = os.path.basename(path[:-14])
    print name

    out = re.search(r'at[0-9]+-[0-9,I, N, F]+', name)
    label = out.group(0)[2:]
    print label

    tree[label] = root2array(path)
    # branches:
    # ['Nev','tof_reco', 'PID', 'tof_gen',
    # 'P_reco', 'vtx_SumPT2', 'vtx_NDOF','vtx_SumPT',
    # 'Zout','Tout','pt','ctgtheta','phi','d0','dz', 'L',
    # 'sigma_pt', 'sigma_d0', 'sigma_dz', 'sigma_Tin',
    # 'M_reco', 'beta_reco'
    # ]
    arr_pt_vtxpt = np.column_stack((tree[label]['vtx_SumPT'], tree[label]['pt']))
    h_pt[label] = create_TH2D(arr_pt_vtxpt, 'h_pt_'+name, 'SM particle generated tracks, HardQCD #hat{p}_{T} '+label+' GeV', binning=binning, axis_title=axis_title)



c_pt = rt.TCanvas('c_'+name, 'c_'+name, 1600, 600)
c_pt.Divide(len(h_pt.keys())+1,1)

line = rt.TLine()
line.SetLineColor(2)
line.SetLineWidth(2)
line.SetLineStyle(7)

for i, k in enumerate(h_pt.keys()):
    c_pt.cd(i+1)
    h_pt[k].Draw('colz')
    line.DrawLine(min_vtx_pt, min_trk_pt, min_vtx_pt, binning[5])
    line.DrawLine(min_vtx_pt, min_trk_pt, binning[2], min_trk_pt)

c_pt.Update()
#
#
binning = [40,0,400]
axis_title = ['M_{reco} [GeV]', '# Events / {:.0f} GeV [a.u.]'.format((binning[2]-binning[1])/binning[0])]
h_mass = {}
for i, k in enumerate(h_pt.keys()):
    kin_sel = np.logical_and(np.abs(tree[k]['pt']) >= min_trk_pt, np.abs(tree[k]['vtx_SumPT']) >= min_vtx_pt)
    aux = tree[k]['M_reco'][kin_sel]
    h_mass[k] = create_TH1D(aux, 'h_mass_'+k, k, binning=binning, axis_title=axis_title)
    h_mass[k].SetLineColor(2+2*i)
    h_mass[k].SetLineWidth(2)
#
#     aux = t['M_reco'][np.logical_and(np.logical_not(bsm_sel), kin_sel)]
#     h_mass_sm = create_TH1D(aux, 'h_m_sm_'+name, 'SM',binning=binning, axis_title=axis_title)
#     h_mass_sm.SetLineColor(4)
#     h_mass_sm.SetLineWidth(2)
#
# h_m_stack = rt.THStack('h_stack_m_'+name, 'Stack distribution of reco Mass')
#     if h_mass_sm.GetEntries()>h_mass_bsm.GetEntries():
#         h_m_stack.Add(h_mass_bsm)
#         h_m_stack.Add(h_mass_sm)
#     else:
#         h_m_stack.Add(h_mass_sm)
#         h_m_stack.Add(h_mass_bsm)
#
#     leg_m = rt.TLegend(0.75,0.7, 0.95, 0.95)
#     leg_m.AddEntry(h_mass_sm,'SM','le')
#     leg_m.AddEntry(h_mass_bsm,'BSM','le')
#
c_pt.cd(3)
h_mass['200-INF'].Draw('E1L')
#     h_m_stack.Draw('E1L')
#     h_m_stack.GetHistogram().SetXTitle(axis_title[0])
#     h_m_stack.GetXaxis().SetNdivisions(507)
#     h_m_stack.GetHistogram().SetYTitle(axis_title[1])
#     leg_m.Draw()
rt.gPad.SetLogy()

# ---- Fit the mass distribution -----------
h_sum = rt.TH1D(h_mass['200-INF'])
# h_sum.Add(h_mass_sm, h_mass_bsm)
low_b = 50

'''
[0] = Norm
[1] = a1
[2] = a2
[3] = mu
[4] = sigma
'''
t_str = '((x-[3])/[4])'

left_fun = '(TMath::Sign(0.5,-{}-[1])+0.5)*exp(0.5*[1]*[1])*exp([1]*{})'.format(t_str,t_str)
right_fun = "(TMath::Sign(0.5,{}-[2])+0.5)*exp(0.5*[2]*[2])*exp(-[2]*{})".format(t_str, t_str)
central_fun = '(TMath::Sign(0.5,{}+[1])+0.5)*(TMath::Sign(0.5,[2]-{})+0.5)*exp(-0.5*{}*{})'.format(t_str, t_str, t_str, t_str)


f_fit = rt.TF1("f_fit",'[0]*(' + left_fun + " + " + central_fun + ' + ' + right_fun + ')', low_b , binning[2])
f_fit.Identifier = "GausDoubleExp"

integral = h_sum.GetEntries()*h_sum.GetBinWidth(1)
integral /= np.sqrt(2*np.pi)*h_sum.GetRMS()
f_fit.SetParameter(0, integral) #Norm
f_fit.SetParLimits(0,0.1,1e8)

f_fit.SetParameter(1,0.7) #left alpha
f_fit.SetParLimits(1,0.1,5)

f_fit.SetParameter(2,0.7) #right alpha
f_fit.SetParLimits(2,0.1,5)

f_fit.SetParameter(3,h_sum.GetMean()) #mean
f_fit.SetParLimits(3,binning[1], binning[2])

f_fit.SetParameter(4, h_sum.GetRMS()) #sigma
f_fit.SetParLimits(4, 10, 1000)


f_fit = rt.TF1("f_fit",'expo', low_b , binning[2])
f_fit.Identifier = "Exp"

fit_results = h_sum.Fit(f_fit, 'IL0SR')
chi2_free = [f_fit.GetChisquare(), f_fit.GetNDF()]
print chi2_free

f_fit.SetLineColor(8)
f_fit.SetLineStyle(7)
f_fit.SetLineWidth(2)
f_fit.Draw('same')
#
#     pars.append((stopmass,
#                 f_fit.GetParameter(1),
#                 f_fit.GetParError(1),
#                 f_fit.GetParameter(2),
#                 f_fit.GetParError(2),
#                 f_fit.GetParameter(3),
#                 f_fit.GetParError(3),
#                 f_fit.GetParameter(4),
#                 f_fit.GetParError(4),
#                 f_fit.GetChisquare(),
#                 f_fit.GetNDF()
#                 ))
#
#
c_pt.Update()

f = rt.TFile('_root/histos/RecoMassSpectrum_bkg.root', 'RECREATE')
f.cd()
h_sum.Write()
f.Close()
#     c_pt.SaveAs('~/Desktop/SelectionCutsAndModeling'+name+'.png')
#     donotdelete.append(c_pt)
#
#
# pars_names = ['M_st', 'a1', 'a1_err', 'a2', 'a2_err', 'mu', 'mu_err', 'sigma', 'sigma_err', 'chi2', 'NDF']
# pars_labels = [r'$\alpha_{L}$', r'$\alpha_{R}$', r'$\mu$ [GeV]', r'$\sigma$ [GeV]']
# pdeg = [1,1,1,1]
# N_par = (len(pars_names)-3)/2
# pars = np.array(pars, dtype=zip(pars_names,['<f8']*len(pars_names)))
#
# fig, axarr = plt.subplots(N_par + 1, sharex=True, figsize=(6,12))
# fitted_pars = {}
# for i in range(N_par):
#     pn = pars_names[1+2*i]
#     pen = pars_names[2+2*i]
#     axarr[i].errorbar(pars['M_st'], pars[pn], yerr=pars[pen], fmt='.')
#     axarr[i].set_ylabel(pars_labels[i])
#     axarr[i].grid()
#
#     fitted_pars[pn] = np.polyfit(pars['M_st'], pars[pn], pdeg[i], w=1/np.square(pars[pen]))
#     res = np.polyval(fitted_pars[pn], pars['M_st']) - pars[pn]
#     chi2 = np.sum(np.square(res)/np.square(pars[pen]))
#     # print pn, ':', chi2, ' - ', len(pars['M_st'])-pdeg[i]
#
#     x = np.arange(100, 2200)
#     y = np.polyval(fitted_pars[pn], x)
#     axarr[i].plot(x,y,'--')
#
#     axarr[i].text(np.min(pars['M_st']), np.max(pars[pn]), r'Deg: {}, $\chi^{{2}}={:.1f}$ ({}) '.format(pdeg[i], chi2, len(pars['M_st'])-pdeg[i]))
#
#
#
# axarr[N_par].errorbar(pars['M_st'], pars['NDF'], yerr=np.sqrt(pars['NDF']), fmt='.')
# axarr[N_par].plot(pars['M_st'], pars['chi2'], 'o')
# axarr[N_par].set_ylabel(r'$\chi^{2}$')
# axarr[N_par].set_xlabel(r'$M_{stop}$ [GeV]')
# axarr[N_par].grid()
# axarr[0].set_title('Fit results')
#
# fig.show()
# fig.savefig(os.environ['HOME']+'/Desktop/SignalParametrization.png')

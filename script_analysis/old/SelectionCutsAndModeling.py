import numpy as np
import ROOT as rt
from root_numpy import root2array
from histo_utilities import create_TH2D, create_TH1D
from cebefo_style import cebefo_style
import os, sys, re
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
# plt.style.use('olmo')

cebefo_style()
c_light = 2.99792458E8 #m/s
donotdelete = []

min_trk_pt = 100
min_vtx_pt = 370

file_path = list(sys.argv[1:])
print len(sys.argv)-1

pars = []
xsec_list = []
tracks_per_evt = []

for path in file_path:
    name = os.path.basename(path[:-14])
    print name

    out = re.search('_M[0-9]+_', name)
    out = out.group(0)
    stopmass = int(out[2:-1])

    # if stopmass >300 : continue

    t = root2array(path)
    # branches:
    # ['Nev','tof_reco', 'PID', 'tof_gen',
    # 'P_reco', 'vtx_SumPT2', 'vtx_NDOF','vtx_SumPT',
    # 'Zout','Tout','pt','ctgtheta','phi','d0','dz', 'L',
    # 'sigma_pt', 'sigma_d0', 'sigma_dz', 'sigma_Tin',
    # 'M_reco', 'beta_reco'
    # ]

    c_pt = rt.TCanvas('c_'+name, 'c_'+name, 1600, 600)
    c_pt.Divide(3,1)

    bsm_sel = np.logical_and(np.abs(t['PID']) >= 1000612, np.abs(t['PID']) <= 1093334)
    kin_sel = np.logical_and(np.abs(t['pt']) >= min_trk_pt, np.abs(t['vtx_SumPT']) >= min_vtx_pt)

    arr_pt_vtxpt = np.column_stack((t['vtx_SumPT'], t['pt']))
    axis_title = ['reco vtx sum p_{T} [GeV]', 'reco p_{T} [GeV]']
    binning = [100, 0, 3000, 100, 0, 1500]

    h_BSM = create_TH2D(arr_pt_vtxpt[bsm_sel], 'h_pt_BSM_'+name, 'BSM particle generated tracks, #tilde{{t}}_{{1}} mass {} GeV'.format(stopmass), binning=binning, axis_title=axis_title)
    h_SM = create_TH2D(arr_pt_vtxpt[np.logical_not(bsm_sel)], 'h_pt_SM_'+name, 'SM particle generated tracks', binning=binning, axis_title=axis_title)

    line = rt.TLine()
    line.SetLineColor(2)
    line.SetLineWidth(2)
    line.SetLineStyle(7)


    c_pt.cd(1)
    h_BSM.Draw('colz')
    line.DrawLine(min_vtx_pt, min_trk_pt, min_vtx_pt, binning[5])
    line.DrawLine(min_vtx_pt, min_trk_pt, binning[2], min_trk_pt)

    c_pt.cd(2)
    h_SM.Draw('colz')
    line.DrawLine(min_vtx_pt, min_trk_pt, min_vtx_pt, binning[5])
    line.DrawLine(min_vtx_pt, min_trk_pt, binning[2], min_trk_pt)



    binning = [130,0,stopmass*2.2]
    axis_title = ['M_{reco} [GeV]', '# Events / {:.0f} GeV [a.u.]'.format((binning[2]-binning[1])/binning[0])]

    aux = t['M_reco'][np.logical_and(bsm_sel, kin_sel)]
    h_mass_bsm = create_TH1D(aux, 'h_m_bsm_'+name, 'BSM', binning=binning, axis_title=axis_title)
    h_mass_bsm.SetLineColor(2)
    h_mass_bsm.SetLineWidth(2)

    aux = t['M_reco'][np.logical_and(np.logical_not(bsm_sel), kin_sel)]
    h_mass_sm = create_TH1D(aux, 'h_m_sm_'+name, 'SM',binning=binning, axis_title=axis_title)
    h_mass_sm.SetLineColor(4)
    h_mass_sm.SetLineWidth(2)

    h_m_stack = rt.THStack('h_stack_m_'+name, 'Stack distribution of reco Mass')
    if h_mass_sm.GetEntries()>h_mass_bsm.GetEntries():
        h_m_stack.Add(h_mass_bsm)
        h_m_stack.Add(h_mass_sm)
    else:
        h_m_stack.Add(h_mass_sm)
        h_m_stack.Add(h_mass_bsm)

    leg_m = rt.TLegend(0.75,0.7, 0.95, 0.95)
    leg_m.AddEntry(h_mass_sm,'SM','le')
    leg_m.AddEntry(h_mass_bsm,'BSM','le')

    c_pt.cd(3)
    h_m_stack.Draw('E1L')
    h_m_stack.GetHistogram().SetXTitle(axis_title[0])
    h_m_stack.GetXaxis().SetNdivisions(507)
    h_m_stack.GetHistogram().SetYTitle(axis_title[1])
    leg_m.Draw()
    rt.gPad.SetLogy()

    # ---- Fit the mass distribution -----------
    h_sum = rt.TH1D(h_mass_bsm)
    # h_sum.Add(h_mass_sm, h_mass_bsm)

    '''
    [0] = Norm
    [1] = aL
    [2] = aR
    [3] = mu
    [4] = sigma
    '''
    t_str = '((x-[3])/[4])'

    left_fun = '(TMath::Sign(0.5,-{}-[1])+0.5)*exp(0.5*[1]*[1])*exp([1]*{})'.format(t_str,t_str)
    right_fun = "(TMath::Sign(0.5,{}-[2])+0.5)*exp(0.5*[2]*[2])*exp(-[2]*{})".format(t_str, t_str)
    central_fun = '(TMath::Sign(0.5,{}+[1])+0.5)*(TMath::Sign(0.5,[2]-{})+0.5)*exp(-0.5*{}*{})'.format(t_str, t_str, t_str, t_str)

    # low_b = 0
    # if stopmass>750:
    #     low_b = 300
    # else:
    #     low_b = binning[1]

    aux = t['M_reco'][np.logical_and(bsm_sel, kin_sel)]

    f_fit = rt.TF1("f_fit",'[0]*(' + left_fun + " + " + central_fun + ' + ' + right_fun + ')', np.percentile(aux, 3) , np.percentile(aux, 97))
    f_fit.Identifier = "GausDoubleExp"

    integral = h_sum.GetEntries()*h_sum.GetBinWidth(1)
    integral /= np.sqrt(2*np.pi)*h_sum.GetRMS()
    f_fit.SetParameter(0, integral) #Norm
    f_fit.SetParLimits(0,0.1,1e8)

    f_fit.SetParameter(1,1.7) #left alpha
    f_fit.SetParLimits(1,0.1,5)

    f_fit.SetParameter(2,1.7) #right alpha
    f_fit.SetParLimits(2,0.1,5)

    f_fit.SetParameter(3,h_sum.GetMean()) #mean
    f_fit.SetParLimits(3,binning[1], binning[2])

    f_fit.SetParameter(4, h_sum.GetRMS()*0.8) #sigma
    f_fit.SetParLimits(4, 1, 1000)

    f_fit.SetLineColor(8)
    f_fit.SetLineStyle(7)
    f_fit.SetLineWidth(2)

    fit_results = h_sum.Fit(f_fit, 'IL0SQR')
    chi2_free = [fit_results.Chi2(), fit_results.Ndf()]

    f_fit.Draw('same')
    c_pt.Update()
    c_pt.SaveAs('~/Desktop/SelectionCutsAndModeling'+name+'.png')
    donotdelete.append(c_pt)

    # Get the xsec
    xsec_files = os.path.dirname(path) + '/xsec.txt'
    fxsec = open(xsec_files)
    err = []
    xsec = []
    for l in fxsec.readlines():
        out = re.search(r'[0-9]+.[0-9]+e-[0-9]+  [0-9]+.[0-9]+e-[0-9]+', l)
        out =  str(out.group(0))
        out = out.split('  ')
        xsec.append(float(out[0])*1e9)
        err.append(float(out[1])*1e9)

    err2 = np.square(np.array(err))
    xsec = np.array(xsec)
    xsec = np.sum(xsec*err2/np.sum(err2))
    err = 1./np.sqrt(np.sum(1./err2))

    print 'XSec = {:.2g} +/- {:.1g}'.format(xsec, err)
    print 'Tracks = ', h_sum.GetEntries()
    out = re.search(r'_[0-9]+k', name)
    out = float(out.group(0)[1:-1])*1000
    print 'Events gen = ', out

    xsec_list.append(xsec)
    tracks_per_evt.append(h_sum.GetEntries()/float(out))
    print 'eff = ', tracks_per_evt[-1]

    pars.append((stopmass,
                f_fit.GetParameter(1),
                f_fit.GetParError(1),
                f_fit.GetParameter(2),
                f_fit.GetParError(2),
                f_fit.GetParameter(3),
                f_fit.GetParError(3),
                f_fit.GetParameter(4),
                f_fit.GetParError(4),
                fit_results.Chi2(),
                fit_results.Ndf()
                ))

    print '\n\n'




pars_names = ['M_st', 'aL', 'aL_err', 'aR', 'aR_err', 'mu', 'mu_err', 'sigma', 'sigma_err', 'chi2', 'NDF']
pars_labels = [r'$\alpha_{L}$', r'$\alpha_{R}$', r'$\mu$ [GeV]', r'$\sigma$ [GeV]']
pdeg = [1,1,1,3]
N_par = (len(pars_names)-3)/2
pars = np.array(pars, dtype=zip(pars_names,['<f8']*len(pars_names)))

fig = plt.figure(3)
plt.plot(pars['M_st'], tracks_per_evt, '*')

eff_pred = interp1d(pars['M_st'], tracks_per_evt)
aux = np.arange(np.min(pars['M_st']), np.max(pars['M_st']), 1.)
plt.plot(aux, eff_pred(aux), '--')

plt.xlabel(r'$M_{stop}$ [GeV]')
plt.ylabel(r'$\epsilon = $ tks/evt')

plt.title('Efficiency for signal event')
# plt.yscale("log", nonposy='clip')
# plt.xscale("log")
plt.grid()
fig.savefig(os.environ['HOME']+'/Desktop/Efficiency.png')
print 'Masses:', np.sort(np.array(pars['M_st']))
print 'Eff:', np.sort(tracks_per_evt)

fig = plt.figure(2)
plt.plot(pars['M_st'], xsec_list, '*')
coeff = np.polyfit(pars['M_st'], xsec_list, 2)
aux = np.arange(np.min(pars['M_st']), np.max(pars['M_st']), 1.)
eff_pred = np.polyval(coeff, aux)
plt.plot(aux, eff_pred, '--')
plt.xlabel(r'$M_{stop}$ [GeV]')
plt.ylabel(r'xsec [pb]')
title = 'Fit: '
for i, p in enumerate(coeff):
    title += ' {:1.2e}'.format(p)
plt.title(title)
plt.yscale("log", nonposy='clip')
plt.xscale("log")
fig.savefig(os.environ['HOME']+'/Desktop/XSec.png')

fig = plt.figure(1 )
fig, axarr = plt.subplots(N_par + 1, sharex=True, figsize=(6,12))
fitted_pars = {}
f_pars = open('/Users/olmo/Desktop/fit_pars.txt', 'w')
for i in range(N_par):
    pn = pars_names[1+2*i]
    pen = pars_names[2+2*i]
    axarr[i].errorbar(pars['M_st'], pars[pn], yerr=pars[pen], fmt='.')
    axarr[i].set_ylabel(pars_labels[i])
    axarr[i].grid()

    fitted_pars[pn] = np.polyfit(pars['M_st'], pars[pn], pdeg[i], w=1/np.square(pars[pen]))
    res = np.polyval(fitted_pars[pn], pars['M_st']) - pars[pn]
    chi2 = np.sum(np.square(res)/np.square(pars[pen]))
    # print pn, ':', chi2, ' - ', len(pars['M_st'])-pdeg[i]
    line2write = pn
    for p in fitted_pars[pn]:
        line2write += '  {:1.2e}'.format(p)
    line2write += '\n'
    f_pars.write(line2write)

    x = np.arange(0.9*np.min(pars['M_st']), 1.1*np.max(pars['M_st']))
    y = np.polyval(fitted_pars[pn], x)
    axarr[i].plot(x,y,'--')

    axarr[i].text(np.min(pars['M_st']), np.max(pars[pn]), r'Deg: {}, $\chi^{{2}}={:.1f}$ ({}) '.format(pdeg[i], chi2, len(pars['M_st'])-pdeg[i]))
f_pars.close()


axarr[N_par].errorbar(pars['M_st'], pars['NDF'], yerr=np.sqrt(pars['NDF']), fmt='.')
axarr[N_par].plot(pars['M_st'], pars['chi2'], 'o')
axarr[N_par].set_ylabel(r'$\chi^{2}$')
axarr[N_par].set_xlabel(r'$M_{stop}$ [GeV]')
axarr[N_par].grid()
axarr[0].set_title('Fit results')

fig.show()
fig.savefig(os.environ['HOME']+'/Desktop/SignalParametrization.png')

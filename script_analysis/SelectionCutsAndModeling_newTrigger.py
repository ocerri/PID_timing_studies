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

min_trk_pt = 30
max_dm = 0.1
min_vtx_pt = 100

binning_vtx = [100, 0, 3000, 50, 0.001, 1]
binning_pt = [100, 0, 800, 50, 0.001, 1]


file_path = list(sys.argv[1:])
print len(sys.argv)-1

pars = []
efficiency = []


for path in file_path:
    name = os.path.basename(path[:-14])
    print name

    out = re.search('_M[0-9]+_', name)
    out = out.group(0)
    stopmass = int(out[2:-1])
    if stopmass > 1000 or stopmass < 100:
        continue

    t = root2array(path)
    dm = t['mh']/t['ml'] - 1

    arr_dm_vtxpt = np.column_stack((t['vtx_SumPT'][t['ml']>0], dm[t['ml']>0]))
    h_vtx = create_TH2D(arr_dm_vtxpt, 'h_dm_vtxpt_'+name, name, binning=binning_vtx, axis_title=['#sum p_{T} [GeV]','m_{H}/m_{L} - 1'])

    arr_dm_pt = np.column_stack((t['pth'][t['ml']>0], dm[t['ml']>0]))
    h_pt = create_TH2D(arr_dm_pt, 'h_dm_pt_'+name, name, binning=binning_pt, axis_title=['p_{T}^{H} [GeV]','m_{H}/m_{L} - 1'])


    kin_sel = np.logical_and(np.abs(t['pth']) >= min_trk_pt, np.abs(t['vtx_SumPT']) >= min_vtx_pt)
    kin_sel = np.logical_and(kin_sel, t['ml']>0)
    kin_sel = np.logical_and(kin_sel, dm<max_dm)

    binning_mass = [50,stopmass/2.,stopmass*2.]
    arr_mh = (t['mh'][kin_sel] + t['ml'][kin_sel])/2.
    h_mass = create_TH1D(arr_mh, 'h_mass_'+name, name, binning=binning_mass, axis_title=['(m_{{H}} + m_{{L}})/2 [GeV]','Events / {:.1} GeV'.format(float(binning_mass[2]-binning_mass[1])/binning_mass[0])])
    # ---- Fit the mass distribution -----------
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

    f_fit = rt.TF1("f_fit",'[0]*(' + left_fun + " + " + central_fun + ' + ' + right_fun + ')', np.percentile(arr_mh, 0.5) , np.percentile(arr_mh, 99.5))
    f_fit.Identifier = "GausDoubleExp"

    integral = h_mass.GetEntries()*h_mass.GetBinWidth(1)
    integral /= np.sqrt(2*np.pi)*h_mass.GetRMS()
    f_fit.SetParameter(0, integral) #Norm
    f_fit.SetParLimits(0,0.1,1e8)

    f_fit.SetParameter(1,1.7) #left alpha
    f_fit.SetParLimits(1,0.1,5)
    # f_fit.FixParameter(1,5) #left alpha

    f_fit.SetParameter(2,1.7) #right alpha
    f_fit.SetParLimits(2,0.1,5)
    # f_fit.FixParameter(2,5) #right alpha

    f_fit.SetParameter(3,h_mass.GetMean()) #mean
    f_fit.SetParLimits(3,binning_mass[1], binning_mass[2])

    f_fit.SetParameter(4, h_mass.GetRMS()*0.8) #sigma
    f_fit.SetParLimits(4, 1, 1000)

    fit_results = h_mass.Fit(f_fit, 'IL0SQR')
    chi2_free = [fit_results.Chi2(), fit_results.Ndf()]

    print 'Tracks = ', h_mass.GetEntries()
    out = re.search(r'_[0-9]+k', name)
    out = float(out.group(0)[1:-1])*1000
    print 'Events gen = ', out

    efficiency.append(h_mass.GetEntries()/float(out))
    print 'eff = {:.4f}\n\n'.format(efficiency[-1])

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

    '''------------- Drawing -------------------'''
    c_mass = rt.TCanvas('c_mass_res_'+name, 'c_mass_res_'+name, 1600, 600)
    c_mass.Divide(3,1)

    line = rt.TLine()
    line.SetLineColor(2)
    line.SetLineWidth(2)
    line.SetLineStyle(7)

    c_mass.cd(1)
    h_vtx.DrawCopy('colz')
    # rt.gPad.SetLogy()
    line.DrawLine(min_vtx_pt, binning_vtx[4], min_vtx_pt, max_dm)
    line.DrawLine(min_vtx_pt, max_dm, binning_vtx[2], max_dm)

    c_mass.cd(2)
    h_pt.DrawCopy('colz')
    # rt.gPad.SetLogy()
    line.DrawLine(min_trk_pt, binning_pt[4], min_trk_pt, max_dm)
    line.DrawLine(min_trk_pt, max_dm, binning_pt[2], max_dm)

    c_mass.cd(3)
    rt.gPad.SetLogy()
    h_mass.DrawCopy('E1')
    f_fit.SetLineColor(2)
    f_fit.SetLineStyle(7)
    f_fit.SetLineWidth(3)
    f_fit.Draw('same')

    c_mass.Update()
    donotdelete.append((c_mass))
    # c_mass.SaveAs('/Users/olmo/Desktop/SelectionCutsAndModeling_2cat'+name+'.png')


pars_names = ['M_st', 'aL', 'aL_err', 'aR', 'aR_err', 'mu', 'mu_err', 'sigma', 'sigma_err', 'chi2', 'NDF']
pars_labels = [r'$\alpha_{L}$', r'$\alpha_{R}$', r'$\mu$ [GeV]', r'$\sigma$ [GeV]']
pdeg = [1,1,1,3]
N_par = (len(pars_names)-3)/2
pars = np.array(pars, dtype=zip(pars_names,['<f8']*len(pars_names)))


print 'Masses:', '{'+', '.join(['{}'.format(int(e)) for e in np.sort(np.array(pars['M_st']))]) + '}'
print 'Eff:', '{'+', '.join(['{:1.2e}'.format(e) for e in np.sort(efficiency)]) + '}'


fig = plt.figure()
fig, axarr = plt.subplots(N_par + 2, sharex=True, figsize=(6,12))
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
axarr[N_par].grid()
# axarr[0].set_title('Fit results')

axarr[N_par+1].plot(np.sort(pars['M_st']), np.sort(efficiency), '*--')
axarr[N_par+1].set_xlabel(r'$M_{stop}$ [GeV]')
axarr[N_par+1].set_ylabel(r'$\epsilon$')
# axarr[N_par+1].set_yscale("log")
axarr[N_par+1].grid()

fig.show()
fig.savefig(os.environ['HOME']+'/Desktop/SignalParametrization.png')

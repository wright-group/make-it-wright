import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import pearsonr
import WrightTools as wt
import makeitwright.spectra as spectra, styles
from .helpers import norm, roi


pi = np.pi


def fromBruker(*filepaths):
    d = []
    for filepath in filepaths:
        dtype = "Locked Coupled"
        header_size=None
        with open(filepath) as f:
            txt = f.readlines()
        for i, line in enumerate(txt):
            if "ScanType" in line:
                dtype = line.split('=')[-1].strip()
            if "[Data]" in line:
                header_size = i+2
        if header_size is None:
            try:
                arr = np.genfromtxt(filepath, skip_header=166, delimiter=',')
                print("Data header was not identified in file. Data in instance may not reflect complete file information.")
            except:
                print("Unable to read data from file due to lack of expected data header.")
        else:
            arr = np.genfromtxt(filepath, skip_header=header_size, delimiter=',')
        
        if arr.size > 0:
            deg_arr = arr[:,0].flatten()
            ch_arr = arr[:,1].flatten()
            pat = wt.Data(name=filepath.split('/')[-1])
            pat.create_channel('sig', values=ch_arr)
            pat.create_channel('norm', values=norm(ch_arr, 1, 100))
            pat.create_channel('log', values=np.log(norm(ch_arr, 1, 100)))
            if dtype=="Locked Coupled":
                pat.create_variable('ang', values=deg_arr, units='deg')
                pat.transform('ang')
                pat.attrs['acquisition'] = 'XRD_2theta'
            if dtype=="Z-Drive":
                pat.create_variable('z', values=deg_arr, units='mm')
                pat.transform('z')
                pat.attrs['acquisition'] = 'XRD_2theta'
            pat.attrs['dtype'] = 'spectrum'
            d.append(pat)
        else:
            print(f'file {filepath} was loaded but had no values')
            
    return d

def get_fits(data, channel='norm', function='gauss', xrange='all'):
    def gauss(x, a, u, s):
        return a*np.exp(-((x-u)/(2*s))**2)
    def cauchy(x, a, u, s):
        return a/(pi*s*(1+((x-u)/s)**2))     
    
    functions = {
        'gauss' : gauss,
        'cauchy' : cauchy,
        'lorentz' : cauchy
        }
    
    fits = {}
    for i in range(len(data)):
        if xrange != 'all':
            out = roi(data[i], {'ang':xrange})
        else:
            out = data[i]
        ch_arr = norm(out[channel][:],0,1)
        fit, cov  = curve_fit(functions[function], out['ang'][:], ch_arr, bounds=([0.01,xrange[0],0.001],[1,xrange[1],0.03]), maxfev=1000*len(out['ang'][:]))
        std = np.sqrt(np.diag(cov))
            
        fitdict = {
            'A' : (fit[0], std[0]),
            'u' : (fit[1], std[1]),
            's' : (fit[2], std[2])
            }
    
        fitdict['r2'] = pearsonr(ch_arr, functions[function](out['ang'][:], *fit))[0]**2
            
        fits[f'{i}_{data[i].natural_name}'] = fitdict
        
        fig, gs = wt.artists.create_figure()
        ax = plt.subplot(gs[0])
        ax.plot(out['ang'][:], functions[function](out['ang'][:], *fit), linewidth=3, color='red')
        ax.scatter(out['ang'][:], ch_arr, color='black', s=10)
        ax.set_title(f'{data[i].natural_name[-20:]}')
        ax.set_ylabel("intensity (a.u.)")
        ax.set_xlabel("diffraction angle (deg. 2\u03B8)")
        ax.set_ylim(0,1)

    return fits

def plot_patterns(data, **kwargs):
    params = {}
    params.update(styles.spectra_XRD_pattern)
    params.update(**kwargs)
    spectra.plot_spectra(data, **params)
    
def fig1(c1, m1, s1, e1, c2, m2, s2, e2, figsize=(6.7,1.35), xrange=[0,100], vrange=(3.75,4.05), colors=['red','blue'], linewidth=0.5, elinewidth=1, reference_lines=[10,20,30,40,50,60,70,80,90]):
    #setup plot frame
    fig, ax = plt.subplots(figsize=figsize)
    
    #plot data
    ax.scatter(c1,m1,s=s1,linewidth=linewidth,color=colors[0], alpha=1)
    ax.scatter(c2,m2,s=s2,linewidth=linewidth,color=colors[1], alpha=1)
    ax.errorbar(c1, m1, yerr=e1, color='grey', linewidth=linewidth, linestyle=':', elinewidth=elinewidth, capsize=2, zorder=0)
    ax.errorbar(c2, m2, yerr=e2, color='grey', linewidth=linewidth, linestyle=':', elinewidth=elinewidth, capsize=2, zorder=0)

    if reference_lines is not None:
        if type(reference_lines) is not list:
            reference_lines = [reference_lines]
        for line in reference_lines:
            ax.axvline(x=line, zorder=0, linewidth=0.25, color='grey', alpha=0.5)
    
    #adjust plot frame
    ax.set_xlim(*xrange)
    ax.set_ylim(*vrange)
    ax.set_xlabel("")
    ax.set_ylabel("")
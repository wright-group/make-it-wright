import numpy as np
import matplotlib.pyplot as plt
import WrightTools as wt
from scipy.optimize import curve_fit
from scipy.stats import pearsonr

from .core.helpers import get_axes, get_channels, set_label, roi
from .core import spectra
from .core import styles




def get_fits(data, channel='norm', function='biexp'):
    def exp(t, a, td):
        return a*np.exp(-t/td)     
    def biexp(t, a1, td1, a2, td2):
        return a1*np.exp(-t/td1) + a2*np.exp(-t/td2)
    
    functions = {
        "exp" : exp,
        "biexp" : biexp
        }
    
    fits = {}
    for i in range(len(data)):
        out = roi(data[i], {'t':[0]})
        fit, cov  = curve_fit(functions[function], out['t'][:], out[channel][:], bounds=(0,1000000), maxfev=1000*len(out['t'][:]))
        std = np.sqrt(np.diag(cov))
        
        if function=='biexp':
            sort = sorted(zip([fit[1],fit[3]], [std[1],std[3]], [fit[0],fit[2]], [std[0],std[2]]))
            tsorted = [t for t, tstd, a, astd in sort]
            tstdsorted = [tstd for t, tstd, a, astd in sort]
            asorted = [a for t, tstd, a, astd in sort]
            astdsorted = [astd for t, tstd, a, astd in sort]
            fitdict = {
                'function' : function,
                'A1' : (asorted[0], astdsorted[0]),
                't1' : (tsorted[0], tstdsorted[0]),
                'A2' : (asorted[1], astdsorted[1]),
                't2' : (tsorted[1], tstdsorted[1])
                }
            
        if function=='exp':
            fitdict = {
                'function' : function,
                'A' : (fit[0], std[0]),
                't' : (fit[1], std[1])
                }
        
        fitdict['r2'] = pearsonr(out[channel][:], functions[function](out['t'][:], *fit))[0]**2
            
        fits[f'{i}_{data[i].natural_name}'] = fitdict
        
        fig, gs = wt.artists.create_figure()
        ax = plt.subplot(gs[0])
        ax.scatter(out['t'][:], out[channel][:], color='orange', alpha=0.5, s=3)
        ax.plot(out['t'][:], functions[function](out['t'][:], *fit), '--', linewidth=3, color='blue')
        ax.set_title(f'{data[i].natural_name[:20]}')
        ax.set_ylabel("counts")
        ax.set_xlabel("t (ns)")
        ax.set_yscale('log')
        ax.set_ylim(1,np.max(out[channel][:]))

    return fits

def plot_transients(data, **kwargs):
    params = {}
    params.update(styles.spectra_TRPL)
    params.update(**kwargs)
    
    spectra.plot_spectra(data, **params)

def remove_dark_counts(data, channel, time_range, name=None):
    channel = get_channels(data,channel)[0]
    k = np.average(roi(data,{0:time_range},return_arrs=True)[channel])
    if name is None:
        name = f"{channel}_nodark"
    data.create_channel(name, values=data[channel].points-k)
    set_label(data,name,"PL intensity (a.u.)")

def compute_biexponential_fit(data, channel, axis=0, IRF=None, bounds=None):
    def biexp(t,a1,t1,a2,t2,offset):
        return np.heaviside(t-offset,1) * (a1*np.exp(-(t-offset)/t1)+a2*np.exp(-(t-offset)/t2))
    def convolved_biexp(t,a1,t1,a2,t2,offset,airf):
        return np.convolve((np.heaviside(t-offset,1) * (a1*np.exp(-(t-offset)/t1)+a2*np.exp(-(t-offset)/t2))),airf*IRF,mode='same')
    
    axis, channel = get_axes(data, axis)[0], get_channels(data, channel)[0]
    t_arr, sig_arr = data[axis].points, data[channel].points
    t_arr -= t_arr[np.argmax(sig_arr)]
    maxfev = t_arr.size*1000000
    if bounds is None:
        a1bounds, a2bounds = (0.01*np.max(sig_arr),np.max(sig_arr)), (0.01*np.max(sig_arr),np.max(sig_arr))
        t1bounds, t2bounds = (0,np.max(t_arr)/2), (0,np.max(t_arr))
        airfbounds = (0.04,0.06)
        t_range = np.max(t_arr)-np.min(t_arr)
        offset_bounds = (-0.1*t_range, 0.1*t_range)
        bounds = ([a1bounds[0],t1bounds[0],a2bounds[0],t2bounds[0],offset_bounds[0]],[a1bounds[1],t1bounds[1],a2bounds[1],t2bounds[1],offset_bounds[1]])
    print(f'{bounds}')
    if IRF is not None:
        bounds = ([a1bounds[0],t1bounds[0],a2bounds[0],t2bounds[0],offset_bounds[0],airfbounds[0]],[a1bounds[1],t1bounds[1],a2bounds[1],t2bounds[1],offset_bounds[1],airfbounds[1]])
        fit, cov = curve_fit(convolved_biexp, t_arr, sig_arr, bounds=bounds, maxfev=maxfev)
        fit_arr = convolved_biexp(t_arr, *fit)
        fit_type = "simple biexponential"
    else:
        fit, cov = curve_fit(biexp, t_arr, sig_arr, bounds=bounds, maxfev=maxfev)
        fit_arr = biexp(t_arr, *fit)
        fit_type = "biexponential convolved with IRF"
    err = np.sqrt(np.diag(cov))
    print(f"{fit_arr.size}, {sig_arr.size}")

    fit_dict = {
        'a1' : (fit[0],err[0]),
        't1' : (fit[1],err[1]),
        'a2' : (fit[2],err[2]),
        't2' : (fit[3], err[3]),
        'offset' : (fit[-1],err[-1]),
        'r2': pearsonr(fit_arr,sig_arr)[0]**2
    }
    if IRF is not None:
        fit_dict['aIRF'] = (fit[4],err[4])
    print(f"{fit_dict}")

    data.create_channel(f"fit_{channel}", values=fit_arr)
    set_label(data, f"fit_{channel}", "fit magnitude")
    data[f"fit_{channel}"].attrs['fit'] = str(fit_dict)

    print(f"biexponential fit applied to channel {channel} of data {data.natural_name}:")
    print(f"fit type: {fit_type}")
    print(f"fit parameters:")
    for key, value in fit_dict.items():
        if key != 'r2':
            print(f"    {key} : {value[0]} +/- {value[1]}")
    print(f"    r2 : {fit_dict['r2']}")
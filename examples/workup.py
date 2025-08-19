from numpy import inf
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import pearsonr, moyal, cauchy
import pathlib

from scipy.signal import savgol_filter
from matplotlib import pyplot as plt

import makeitwright as mw


parse = mw.parsers.parse
__at = mw.helpers.find_nearest()
roi = mw.helpers.roi
set_label = mw.helpers.set_label
norm = mw.helpers.norm


def lorentz_fit_2(data, channel='darksub', xrange='all', bounds=None, plot=False):

    if bounds is None:
        bounds = (0,inf)

    def cauchy2(x, a1, u1, s1, a2, u2, s2):
        return a1*cauchy.pdf(x, u1, s1)+a2*cauchy.pdf(x, u2, s2)

    fits = {}

    for i in range(len(data)):
        out = data[i]
        ch_arr = out[channel][:]
        fit, cov  = curve_fit(cauchy2, out['wl'][:], ch_arr, maxfev=1000*len(out['wl'][:]), bounds=bounds)
        std = np.sqrt(np.diag(cov))

        fitdict = {
            'a1' : (fit[0], std[0]),
            'u1' : (fit[1], std[1]),
            's1' : (fit[2], std[2]),
            'a2' : (fit[3], std[3]),
            'u2' : (fit[4], std[4]),
            's2' : (fit[5], std[5])
            }

        fitdict['r2'] = pearsonr(ch_arr, cauchy2(out['wl'][:], *fit))[0]**2

        fits[f'{i}_{data[i].natural_name}'] = fitdict
        
        if plot:
            fig, ax = plt.subplots(figsize=(4,4))
            ax.scatter(out['wl'][:], ch_arr, color='lightgray', s=20)
            ax.plot(out['wl'][:], fit[0]*cauchy.pdf(out['wl'][:], *fit[1:3]), linewidth=1.5, color='cornflowerblue')
            ax.plot(out['wl'][:], fit[3]*cauchy.pdf(out['wl'][:], *fit[4:]), linewidth=1.5, color='cornflowerblue')
            ax.plot(out['wl'][:], cauchy2(out['wl'][:], *fit), linewidth=1.5, color='red')
            ax.set_title(f'{data[i].natural_name[-20:]}')
            ax.set_ylabel("PL intensity (a.u.)")
            ax.set_xlabel("wavelength (nm)")

    return fits

def moyal_fit_2(data, channel='darksub', xrange='all', bounds=None, plot=False):
    
    if bounds is None:
        bounds = (0,inf)
    
    def moyal2(x, a1, u1, s1, a2, u2, s2):
       return a1*moyal.pdf(x, u1, s1)+a2*moyal.pdf(x, u2, s2)
   
    fits = {}

    for i in range(len(data)):
        out = data[i]
        ch_arr = out[channel][:]
        fit, cov  = curve_fit(moyal2, out['wl'][:], ch_arr, maxfev=1000*len(out['wl'][:]), bounds=bounds)
        std = np.sqrt(np.diag(cov))

        fitdict = {
            'a1' : (fit[0], std[0]),
            'u1' : (fit[1], std[1]),
            's1' : (fit[2], std[2]),
            'a2' : (fit[3], std[3]),
            'u2' : (fit[4], std[4]),
            's2' : (fit[5], std[5])
            }

        fitdict['r2'] = pearsonr(ch_arr, moyal2(out['wl'][:], *fit))[0]**2
        
        fits[f'{i}_{data[i].natural_name}'] = fitdict
        
        if plot:
            fig, ax = plt.subplots(figsize=(4,4))
            ax.scatter(out['wl'][:], ch_arr, color='lightgray', s=20)
            ax.plot(out['wl'][:], fit[0]*moyal.pdf(out['wl'][:], *fit[1:3]), linewidth=1.5, color='cornflowerblue')
            ax.plot(out['wl'][:], fit[3]*moyal.pdf(out['wl'][:], *fit[4:]), linewidth=1.5, color='cornflowerblue')
            ax.plot(out['wl'][:], moyal2(out['wl'][:], *fit), linewidth=1.5, color='red')
            ax.set_title(f'{data[i].natural_name[-20:]}')
            ax.set_ylabel("PL intensity (a.u.)")
            ax.set_xlabel("wavelength (nm)")

    return fits

def testfit(data, channel='darksub', xrange='all', bounds=None, plot=False):
    
    if bounds is None:
        bounds = (0,inf)
    
    def test2(x, a1, u1, s1, a2, u2, s2):
       return a1*moyal.pdf(x, u1, s1)+a2*cauchy.pdf(x, u2, s2)
   
    fits = {}

    for i in range(len(data)):
        out = data[i]
        ch_arr = out[channel][:]
        fit, cov  = curve_fit(test2, out['wl'][:], ch_arr, maxfev=1000*len(out['wl'][:]), bounds=bounds)
        std = np.sqrt(np.diag(cov))

        fitdict = {
            'a1' : (fit[0], std[0]),
            'u1' : (fit[1], std[1]),
            's1' : (fit[2], std[2]),
            'a2' : (fit[3], std[3]),
            'u2' : (fit[4], std[4]),
            's2' : (fit[5], std[5])
            }

        fitdict['r2'] = pearsonr(ch_arr, test2(out['wl'][:], *fit))[0]**2
        
        fits[f'{i}_{data[i].natural_name}'] = fitdict
        
        if plot:
            fig, ax = plt.subplots(figsize=(4,4))
            ax.scatter(out['wl'][:], ch_arr, color='lightgray', s=20)
            ax.plot(out['wl'][:], fit[0]*moyal.pdf(out['wl'][:], *fit[1:3]), linewidth=1.5, color='cornflowerblue')
            ax.plot(out['wl'][:], fit[3]*cauchy.pdf(out['wl'][:], *fit[4:]), linewidth=1.5, color='cornflowerblue')
            ax.plot(out['wl'][:], test2(out['wl'][:], *fit), linewidth=1.5, color='red')
            ax.set_title(f'{data[i].natural_name[-20:]}')
            ax.set_ylabel("PL intensity (a.u.)")
            ax.set_xlabel("wavelength (nm)")

    return fits


def residual(a, fit):
    return (a-fit)/a*100


base = pathlib.Path().expanduser().resolve() / r'OneDrive/Documents/UW/research/data local/WG-microscope/biexciton-fluence-dependent-PL_20220909'
fn1 = base / "n1BA"
fn2 = base / 'n2BAMA_CRRsample'
fn3 = base / 'n3BAMA_CRRsample/fluence-series'
fn4 = base / 'n4BAMA'
n1raw, n2raw, n3raw, n4raw = parse(fn1), parse(fn2), parse(fn3), parse(fn4)

for d in n1raw:
    wlmax = __at(d['wl'][:],420)
    ymax = __at(d['y'][:],0.1)
    exposure = float(d.attrs['Exposure Time (secs)'])
    constantbkg = np.average(d['sig'][0:wlmax,:])
    wlbkg = np.average(d['sig'][:,0:ymax],axis=1)
    wlbkg = savgol_filter(wlbkg,101,1)
    wlbkg = wlbkg.reshape((wlbkg.size,1))
    wlbkg = np.repeat(wlbkg,d['y'].size,axis=1)
    constantbkgsub = (d['sig'][:]-constantbkg)/exposure
    wlbkgsub = (d['sig'][:]-wlbkg)/exposure
    d.create_channel('constantbkgsub',values=constantbkgsub)
    d.create_channel('wlbkgsub',values=wlbkgsub)
    set_label(d, 'wlbkgsub', "PL intensity (cps)")
    
for d in n2raw:
    wlmax = __at(d['wl'][:],480)
    ymax = __at(d['y'][:],0.5)
    exposure = float(d.attrs['Exposure Time (secs)'])
    constantbkg = np.average(d['sig'][0:wlmax,:])
    wlbkg = np.average(d['sig'][:,0:ymax],axis=1)
    wlbkg = savgol_filter(wlbkg,101,1)
    wlbkg = wlbkg.reshape((wlbkg.size,1))
    wlbkg = np.repeat(wlbkg,d['y'].size,axis=1)
    constantbkgsub = (d['sig'][:]-constantbkg)/exposure
    wlbkgsub = (d['sig'][:]-wlbkg)/exposure
    d.create_channel('constantbkgsub',values=constantbkgsub)
    d.create_channel('wlbkgsub',values=wlbkgsub)
    set_label(d, 'wlbkgsub', "PL intensity (cps)")

for d in n3raw:
    wlmax = __at(d['wl'][:],540)
    ymax = __at(d['y'][:],0.5)
    exposure = float(d.attrs['Exposure Time (secs)'])
    constantbkg = np.average(d['sig'][0:wlmax,:])
    wlbkg = np.average(d['sig'][:,0:ymax],axis=1)
    wlbkg = savgol_filter(wlbkg,101,1)
    wlbkg = wlbkg.reshape((wlbkg.size,1))
    wlbkg = np.repeat(wlbkg,d['y'].size,axis=1)
    constantbkgsub = (d['sig'][:]-constantbkg)/exposure
    wlbkgsub = (d['sig'][:]-wlbkg)/exposure
    d.create_channel('constantbkgsub',values=constantbkgsub)
    d.create_channel('wlbkgsub',values=wlbkgsub)
    set_label(d, 'wlbkgsub', "PL intensity (cps)")

for d in n4raw:
    wlmax = __at(d['wl'][:],550)
    ymax = __at(d['y'][:],0.5)
    exposure = float(d.attrs['Exposure Time (secs)'])
    constantbkg = np.average(d['sig'][0:wlmax,:])
    wlbkg = np.average(d['sig'][:,0:ymax],axis=1)
    wlbkg = savgol_filter(wlbkg,101,1)
    wlbkg = wlbkg.reshape((wlbkg.size,1))
    wlbkg = np.repeat(wlbkg,d['y'].size,axis=1)
    constantbkgsub = (d['sig'][:]-constantbkg)/exposure
    wlbkgsub = (d['sig'][:]-wlbkg)/exposure
    d.create_channel('constantbkgsub',values=constantbkgsub)
    d.create_channel('wlbkgsub',values=wlbkgsub)
    # set_label(d, 'wlbkgsub', "PL intensity (cps)")

bn1 = roi(n1raw,{'y':'sum'})
bn2 = roi(n2raw,{'y':'sum'})
bn3 = roi(n3raw,{'y':'sum'})
bn4 = roi(n4raw,{'y':'sum'})
yn1 = roi(n1raw,{'wl':'sum'})
yn2 = roi(n2raw,{'wl':'sum'})
yn3 = roi(n3raw,{'wl':'sum'})
yn4 = roi(n4raw,{'wl':'sum'})

for d in bn1+bn2+bn3+bn4:
    a = d['wlbkgsub'][:]
    a = savgol_filter(a,31,1)
    a = norm(a,1,100)
    d.create_channel('norm',values=a)
    set_label(d, 'norm', "norm. PL intensity")
    
for d in yn1+yn2+yn3+yn4:
    a = d['wlbkgsub'][:]
    a = savgol_filter(a,11,1)
    a = norm(a,1,100)
    d.create_channel('norm',values=a)
    set_label(d, 'norm', "norm. PL intensity")

roin1 = roi(bn1,{'wl':[490,610]})
roin2 = roi(bn2,{'wl':[545,665]})
roin3 = roi(bn3,{'wl':[585,705]})
roin4 = roi(bn4,{'wl':[610,730]})

boundsn1 = ([100,515,3,10,528,3],[10000000,527,100,10000000,540,100])
n1lorentz = lorentz_fit_2(roin1,channel='wlbkgsub',xrange='all',bounds=boundsn1)

boundsn2 = ([100,578,3,10,587,3],[10000000,582,100,10000000,597,100])
n2lorentz = lorentz_fit_2(roin2,channel='wlbkgsub',xrange='all',bounds=boundsn2)

boundsn3 = ([500,615,3,10,623,3],[100000000,622,20,100000000,650,20])
n3lorentz = lorentz_fit_2(roin3,channel='wlbkgsub',xrange='all',bounds=boundsn3)

boundsn4 = ([500,646,3,500,655,3],[1000000000,650,20,1000000000,665,20])
n4lorentz = lorentz_fit_2(roin4,channel='wlbkgsub',xrange='all',bounds=boundsn4)

fitsn1 = [n1lorentz[str(i)+'_'+d.natural_name] for i, d in enumerate(bn1)]
fitsn2 = [n2lorentz[str(i)+'_'+d.natural_name] for i, d in enumerate(bn2)]
fitsn3 = [n3lorentz[str(i)+'_'+d.natural_name] for i, d in enumerate(bn3)]
fitsn4 = [n4lorentz[str(i)+'_'+d.natural_name] for i, d in enumerate(bn4)]

pn1uW = np.asarray([float(d.natural_name.split('/')[-1].split('.')[0].split('-')[0]) for d in bn1])+np.asarray([float('0.'+d.natural_name.split('/')[-1].split('.')[0].split('-')[-1]) for d in bn1])
pn2uW = np.asarray([float(d.natural_name.split('/')[-1].split('.')[0].split('-')[0]) for d in bn2])+np.asarray([float('0.'+d.natural_name.split('/')[-1].split('.')[0].split('-')[-1]) for d in bn2])
pn3uW = np.asarray([float(d.natural_name.split('/')[-1].split('.')[0].split('-')[0]) for d in bn3])+np.asarray([float('0.'+d.natural_name.split('/')[-1].split('.')[0].split('-')[-1]) for d in bn3])
pn4uW = np.asarray([float(d.natural_name.split('/')[-1].split('.')[0].split('-')[0]) for d in bn4])+np.asarray([float('0.'+d.natural_name.split('/')[-1].split('.')[0].split('-')[-1]) for d in bn4])
An1, An2, An3, An4 = 1.75E-08, 4.91E-08, 2.85E-08, 8.35E-08
Fn1 = pn1uW*1000/An1/80000000
Fn2 = pn2uW*1000/An2/80000000
Fn3 = pn3uW*1000/An3/80000000
Fn4 = pn4uW*1000/An4/80000000
aXn1 = np.asarray([d['a1'][0] for d in fitsn1])
aXn2 = np.asarray([d['a1'][0] for d in fitsn2])
aXn3 = np.asarray([d['a1'][0] for d in fitsn3])
aXn4 = np.asarray([d['a1'][0] for d in fitsn4])
aXXn1 = np.asarray([d['a2'][0] for d in fitsn1])
aXXn2 = np.asarray([d['a2'][0] for d in fitsn2])
aXXn3 = np.asarray([d['a2'][0] for d in fitsn3])
aXXn4 = np.asarray([d['a2'][0] for d in fitsn4])
uXn1 = np.asarray([d['u1'][0] for d in fitsn1])
uXn2 = np.asarray([d['u1'][0] for d in fitsn2])
uXn3 = np.asarray([d['u1'][0] for d in fitsn3])
uXn4 = np.asarray([d['u1'][0] for d in fitsn4])
uXXn1 = np.asarray([d['u2'][0] for d in fitsn1])
uXXn2 = np.asarray([d['u2'][0] for d in fitsn2])
uXXn3 = np.asarray([d['u2'][0] for d in fitsn3])
uXXn4 = np.asarray([d['u2'][0] for d in fitsn4])
sXn1 = np.asarray([d['s1'][0] for d in fitsn1])
sXn2 = np.asarray([d['s1'][0] for d in fitsn2])
sXn3 = np.asarray([d['s1'][0] for d in fitsn3])
sXn4 = np.asarray([d['s1'][0] for d in fitsn4])
sXXn1 = np.asarray([d['s2'][0] for d in fitsn1])
sXXn2 = np.asarray([d['s2'][0] for d in fitsn2])
sXXn3 = np.asarray([d['s2'][0] for d in fitsn3])
sXXn4 = np.asarray([d['s2'][0] for d in fitsn4])

AXn1, AXXn1 = [], []
for i, d in enumerate(bn1):
    Xfits, XXfits = [uXn1[i], sXn1[i]], [uXXn1[i], sXXn1[i]]
    X = aXn1[i]*cauchy.pdf(d['wl'][:],*Xfits)
    AXn1.append(np.sum(X))
    XX = aXXn1[i]*cauchy.pdf(d['wl'][:],*XXfits)
    AXXn1.append(np.sum(XX))
    d.create_channel('fit', values=X+XX)
    d.create_channel('res', values=residual(d['wlbkgsub'][:],X+XX))
    set_label(d, 'fit', "fit magnitude")
    set_label(d, 'res', "residual (%)")
   
AXn2, AXXn2 = [], []
for i, d in enumerate(bn2):
    Xfits, XXfits = [uXn2[i], sXn2[i]], [uXXn2[i], sXXn2[i]]
    X = aXn2[i]*cauchy.pdf(d['wl'][:],*Xfits)
    AXn2.append(np.sum(X))
    XX = aXXn2[i]*cauchy.pdf(d['wl'][:],*XXfits)
    AXXn2.append(np.sum(XX))
    d.create_channel('fit', values=X+XX)
    d.create_channel('res', values=residual(d['wlbkgsub'][:],X+XX)) 
    set_label(d, 'fit', "fit magnitude")
    set_label(d, 'res', "residual (%)")    

AXn3, AXXn3 = [], []
for i, d in enumerate(bn3):
    Xfits, XXfits = [uXn3[i], sXn3[i]], [uXXn3[i], sXXn3[i]]
    X = aXn3[i]*cauchy.pdf(d['wl'][:],*Xfits)
    AXn3.append(np.sum(X))
    XX = aXXn3[i]*cauchy.pdf(d['wl'][:],*XXfits)
    AXXn3.append(np.sum(XX))
    d.create_channel('fit', values=X+XX)
    d.create_channel('res', values=residual(d['wlbkgsub'][:],X+XX))    
    set_label(d, 'fit', "fit magnitude")
    set_label(d, 'res', "residual (%)")

AXn4, AXXn4 = [], []
for i, d in enumerate(bn4):
    Xfits, XXfits = [uXn4[i], sXn4[i]], [uXXn4[i], sXXn4[i]]
    X = aXn4[i]*cauchy.pdf(d['wl'][:],*Xfits)
    AXn4.append(np.sum(X))
    XX = aXXn4[i]*cauchy.pdf(d['wl'][:],*XXfits)
    AXXn4.append(np.sum(XX))
    d.create_channel('fit', values=X+XX)
    d.create_channel('res', values=residual(d['wlbkgsub'][:],X+XX)) 
    set_label(d, 'fit', "fit magnitude")
    set_label(d, 'res', "residual (%)")

AXn1 = np.asarray(AXn1)
AXXn1 = np.asarray(AXXn1)
AXn2 = np.asarray(AXn2)
AXXn2 = np.asarray(AXXn2)
AXn3 = np.asarray(AXn3)
AXXn3 = np.asarray(AXXn3)
AXn4 = np.asarray(AXn4)
AXXn4 = np.asarray(AXXn4)

aX = [aXn1,aXn2,aXn3,aXn4]
aXX = [aXXn1,aXXn2,aXXn3,aXXn4]
uX = [uXn1,uXn2,uXn3,uXn4]
uXX = [uXXn1,uXXn2,uXXn3,uXXn4]
sX = [sXn1,sXn2,sXn3,sXn4]
sXX = [sXXn1,sXXn2,sXXn3,sXXn4]
AX = [AXn1,AXn2,AXn3,AXn4]
AXX = [AXXn1,AXXn2,AXXn3,AXXn4]

import WrightTools as wt
Xn1 = wt.Data(name='Xn1')
XXn1 = wt.Data(name='XXn1')
Xn2 = wt.Data(name='Xn2')
XXn2 = wt.Data(name='XXn2')
Xn3 = wt.Data(name='Xn3')
XXn3 = wt.Data(name='XXn3')
Xn4 = wt.Data(name='Xn4')
XXn4 = wt.Data(name='XXn4')
Xn1.create_variable('F', values=Fn1)
set_label(Xn1, 'F', "fluence (nJ/cm²)")
XXn1.create_variable('F', values=Fn1)
set_label(XXn1, 'F', "fluence (nJ/cm²)")
Xn2.create_variable('F', values=Fn2)
set_label(Xn2, 'F', "fluence (nJ/cm²)")
XXn2.create_variable('F', values=Fn2)
set_label(XXn2, 'F', "fluence (nJ/cm²)")
Xn3.create_variable('F', values=Fn3)
set_label(Xn3, 'F', "fluence (nJ/cm²)")
XXn3.create_variable('F', values=Fn3)
set_label(XXn3, 'F', "fluence (nJ/cm²)")
Xn4.create_variable('F', values=Fn4)
set_label(Xn4, 'F', "fluence (nJ/cm²)")
XXn4.create_variable('F', values=Fn4)
set_label(XXn4, 'F', "fluence (nJ/cm²)")

X = [Xn1,Xn2,Xn3,Xn4]
XX=[XXn1,XXn2,XXn3,XXn4]

for a, u, s, A, x in zip(aX, uX, sX, AX, X):
    x.create_channel('a', values=a)
    set_label(x, 'a', "magnitude coeff. (a.u.)")
    x.create_channel('u', values=u)
    set_label(x, 'u', "fit median (nm)")
    x.create_channel('s', values=s)
    set_label(x, 's', "scale param. (nm)")
    x.create_channel('A', values=A)
    set_label(x, 'A', "fit area (a.u.)")
    x.transform('F')
    
for a, u, s, A, x in zip(aXX, uXX, sXX, AXX, XX):
    x.create_channel('a', values=a)
    set_label(x, 'a', "magnitude coeff. (a.u.)")
    x.create_channel('u', values=u)
    set_label(x, 'u', "fit median (nm)")
    x.create_channel('s', values=s)
    set_label(x, 's', "scale param. (nm)")
    x.create_channel('A', values=A)
    set_label(x, 'A', "fit area (a.u.)")
    x.transform('F')
    
from scipy.constants import c, h
def Eb(x,xx):
    wlx, wlxx = x['u'][:], xx['u'][:]
    ex, exx = (h*c)/(wlx/1E9)/1.602E-19, (h*c)/(wlxx/1E9)/1.602E-19
    eb = (ex-exx)*1000
    x.create_channel('Eb',values=eb)
    set_label(x, 'Eb', "binding energy (meV)")
Eb(Xn1,XXn1)
Eb(Xn2,XXn2)
Eb(Xn3,XXn3)
Eb(Xn4,XXn4)

for x, xx in zip(X,XX):
    ratio = x['A'][:]/xx['A'][:]
    x.create_channel('ratio',values=ratio)
    
from scipy.stats import linregress, tstd
nrange = [[600,1000],[1000,3000],[5000,10000],[300,700]] #recalibrate these
ratios0 = []
ratios1 = []
ratiosfull = []
for x, xx, n in zip(X,XX,nrange):
    rx0,rxx0 = roi(x,{'F':[0,n[0]]}), roi(xx,{'F':[0,n[0]]})
    rx1,rxx1 = roi(x,{'F':[n[1]]}), roi(xx,{'F':[n[1]]})
    r0 = [linregress(np.log10(rx0['F'][:]),np.log10(rx0['A'][:])),linregress(np.log10(rxx0['F'][:]),np.log10(rxx0['A'][:]))]
    r1 = [linregress(np.log10(rx1['F'][:]),np.log10(rx1['A'][:])),linregress(np.log10(rxx1['F'][:]),np.log10(rxx1['A'][:]))]
    ra = [linregress(np.log10(x['F'][:]),np.log10(x['A'][:])),linregress(np.log10(xx['F'][:]),np.log10(xx['A'][:]))]
    ratios0.append(r0[1][0]/r0[0][0])
    ratios1.append(r1[1][0]/r1[0][0])
    ratiosfull.append(ra[1][0]/ra[0][0])

ebstats = [(np.average(x['Eb'][:]),tstd(x['Eb'][:])) for x in X]

roin1 = roi(bn1,{'wl':[490,610]})
roin2 = roi(bn2,{'wl':[545,665]})
roin3 = roi(bn3,{'wl':[585,705]})
roin4 = roi(bn4,{'wl':[610,730]})
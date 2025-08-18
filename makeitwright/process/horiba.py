import numpy as np
import WrightTools as wt
import makeitwright.styles as styles

from . import spectralprofile
from . import hyperspectral
from . import helpers

def central_wavelength(data):
    pass

def plot_image(data, channel, **kwargs):
    params = {}
    try:
        unit = data['wl'].units
    except KeyError:
        unit = data.constants[0].units
    if unit == 'wn':
        params.update(styles.image_horiba_Raman)
    else:
        params.update(styles.image_horiba_PL)
    params.update(**kwargs)
    
    if len(data.axes) == 3:
        hyperspectral.plot_image(data, channel, **params)
    else:
        spectralprofile.plot_image(data, channel, **params)

def plot_profile(data, channel, profile_axis='y', **kwargs):
    params = {}
    try:
        unit = data['wl'].units
    except KeyError:
        unit = data.constants[0].units

    if data.axes[1].natural_name == 't':
        params.update(styles.profile_horiba_timed_series)
    elif unit == 'wn':
        params.update(styles.profile_horiba_Raman)
    else:
        params.update(styles.profile_horiba_PL)
    params.update(**kwargs)
    
    if len(data.axes) == 3:
        hyperspectral.plot_profile(data, profile_axis, channel, **params)
    else:
        spectralprofile.plot_profile(data, channel, **params)
    
def plot_decomposition(data, channel, **kwargs):
    params = {}
    try:
        unit = data['wl'].units
    except KeyError:
        unit = data.constants[0].units
    if unit == 'wn':
        params.update(styles.decomposition_horiba_Raman)
    else:
        params.update(styles.decomposition_horiba_PL)
    params.update(**kwargs)
    
    if len(data.axes) == 3:
        hyperspectral.plot_decomposition(data, 0, 1, 2, channel, **params)
    else:
        spectralprofile.plot_decomposition(data, 0, 1, channel, **params)

def fromAramis(filepath):
    print("not ready yet, get to work :)")


def horiba_typeID(filepath):
    with open(filepath) as f:
        txt = f.readlines()
    header_size = 0

    for line in txt:
        if "#" in line:
            header_size += 1
            
    wl_arr = np.genfromtxt(filepath, skip_header=header_size, max_rows=1)
    ch_arr = np.genfromtxt(filepath, skip_header=header_size+1)
    xy_cols = ch_arr.shape[1]-wl_arr.shape[0]
    
    dtype = None
    if xy_cols==0:
        dtype = 'LabramHR_spectrum'
    
    if xy_cols==1:
        y = ch_arr[:,0][None,:]
        is_survey = True
        for i in range(1, y.size):
            if y.flatten()[i]-y.flatten()[i-1] != 1:
                is_survey = False
        if is_survey:
            dtype = 'LabramHR_spectrum'
        else:
            dtype = 'LabramHR_linescan'
    
    if xy_cols==2:
        dtype = 'LabramHR_map'
        
    return dtype


def fromLabramHR(filepath, name=None, cps=False):
    if name is None:
        name = filepath.split('/')[-1].split('.')[0]
    with open(filepath) as f:
        txt = f.readlines()
    header_size = 0
    acq = 1
    accum = 1
    spectral_units = 'nm'
    
    for line in txt:
        if "#" in line:
            header_size += 1
        if "Acq. time" in line:
            acq = float(line.split('=\t')[1])
        if "Accumulations" in line:
            accum= int(line.split('=\t')[1])
        if "Range (" in line:
            if 'eV' in line:
                spectral_units='eV'
            elif 'cm' in line:
                spectral_units='wn'
            else:
                spectral_units='nm'
        if 'Spectro' in line:
            if 'eV' in line:
                spectral_units='eV'
            elif 'cm' in line:
                spectral_units='wn'
            else:
                spectral_units='nm'
    total_acq_time = acq*accum
    acq_type = {'wn':"Raman", 'nm':"PL", 'eV':"PL"}
    siglabels = {'wn':"scattering intensity", 'nm':"PL intensity", 'eV':"PL intensity"}
    for key, value in siglabels.items():
        if cps:
            siglabels[key] = siglabels[key] + " (cps)"
        else:
            siglabels[key] = siglabels[key] + " (counts)"
            
    spectlabels = {'wn':"Raman shift (cm\u207b\u2071)", 'nm':"wavelength (nm)", 'eV':"energy (eV)"}

    wl_arr = np.genfromtxt(filepath, skip_header=header_size, max_rows=1)
    ch_arr = np.genfromtxt(filepath, skip_header=header_size+1)
    xy_cols = ch_arr.shape[1]-wl_arr.shape[0]
    
    if xy_cols==0:
        sig = ch_arr[:,1]
        if cps:
            sig = sig/total_acq_time
        wl = ch_arr[:,0]
        d = wt.Data(name=name)
        d.create_variable('wl', values=wl, units=spectral_units)
        d['wl'].label = spectlabels[spectral_units]
        d.create_channel('sig', values=sig)
        d['sig'].label = siglabels[spectral_units]
        d.create_channel('norm', values=helpers.norm(sig, 0, 1))
        d['norm'].label = 'norm. ' + siglabels[spectral_units].split(' (')[0]
        d.transform('wl')
        d.attrs['dtype'] = 'spectrum'
        d.attrs['acquisition'] = 'Horiba_' + acq_type[spectral_units]
        d.attrs['exposure time (s)'] = acq
        d.attrs['number of accumulations'] = accum
        print(f"data from file {filepath.split('/')[-1]} is {acq_type[spectral_units]} spectrum")
        
    if xy_cols==1:
        sig = ch_arr[:,1:].transpose()
        wl = wl_arr[:,None]
        y = ch_arr[:,0][None,:]
        
        is_survey = True
        for i in range(1, y.size):
            if y.flatten()[i]-y.flatten()[i-1] != 1:
                is_survey = False
        
        if is_survey:
            d = []
            for i in range(y.size):
                sig_i = sig[:,i].flatten()
                if cps:
                    sig_i = sig_i/total_acq_time
                spect = wt.Data(name=f"{name}_spect{i}")
                spect.create_variable('wl', values=wl.flatten(), units=spectral_units)
                spect['wl'].label = spectlabels[spectral_units]
                spect.create_channel(name='sig', values=sig_i)
                spect['sig'].label = siglabels[spectral_units]
                spect.create_channel(name='norm', values=helpers.norm(sig_i, 0, 1))
                spect['norm'].label = 'norm. ' + siglabels[spectral_units].split(' (')[0]
                spect.transform('wl')
                spect.attrs['dtype'] = 'spectrum'
                spect.attrs['acquisition'] = 'Horiba_' + acq_type[spectral_units]
                spect.attrs['exposure time (s)'] = acq
                spect.attrs['number of accumulations'] = accum
                d.append(spect)
            print(f"data from file {filepath.split('/')[-1]} is {acq_type[spectral_units]} survey")
        else:
            if cps:
                sig = sig/total_acq_time
            d = wt.Data(name=name)
            d.create_variable('wl', values=wl, units=spectral_units)
            d['wl'].label = spectlabels[spectral_units]
            d.create_channel('sig', values=sig)
            d['sig'].label = siglabels[spectral_units]
            d.create_variable('y', values=y, units='um')
            d['y'].label = "y (µm)"
            d.transform('wl', 'y')
            d.attrs['dtype'] = 'spectralprofile'
            d.attrs['acquisition'] = 'Horiba_' + acq_type[spectral_units]
            d.attrs['exposure time (s)'] = acq
            d.attrs['number of accumulations'] = accum
            print(f"data from file {filepath.split('/')[-1]} is {acq_type[spectral_units]} linescan")
        
    if xy_cols==2:
        xidx = ch_arr[:,0]
        xdim = 1
        for i in range(1,ch_arr.shape[0]):
            if xidx[i] != xidx[i-1]:
                xdim = xdim+1
        ydim = int(ch_arr.shape[0]/xdim)
        
        x  = np.zeros((xdim,1,1))
        y = np.zeros((1,ydim,1))
        wl = wl_arr.reshape([1,1,wl_arr.size])
        sig = np.zeros((xdim,ydim,wl_arr.size))
        
        for i in range(0, ch_arr.shape[0], ydim): x[int(i/ydim),0,0] = ch_arr[i,0]
        y[0,:,0] = ch_arr[:ydim,1]
        for i in range(xdim):
            for j in range(ydim):
                sig[i,j,:] = ch_arr[i*ydim+j,2:].reshape([1,1,wl_arr.size])

        if cps:
            sig = sig/total_acq_time
        d = wt.Data(name=name)
        d.create_channel('sig', values=sig)
        d['sig'].label = siglabels[spectral_units]
        d.create_variable('x', values=x, units='um')
        d['x'].label = "x (µm)"
        d.create_variable('y', values=y, units='um')
        d['y'].label = "y (µm)"
        d.create_variable('wl', values=wl, units=spectral_units)
        d['wl'].label = spectlabels[spectral_units]
        d.transform('x','y','wl')
        d.attrs['dtype'] = 'hyperspectral'
        d.attrs['acquisition'] = 'Horiba_' + acq_type[spectral_units]
        d.attrs['exposure time (s)'] = acq
        d.attrs['number of accumulations'] = accum
        print(f"data from file {filepath.split('/')[-1]} is {acq_type[spectral_units]} map")

    return d

def fromLabramHRTimedSeries(filedir):
    pass
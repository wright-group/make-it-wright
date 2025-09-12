import numpy as np
import WrightTools as wt


def horiba_typeID(filepath):
    """surmise the type of horiba scan (spectrum, linescan, map) without creating data object
    NOTE: dtypes disagree with data attrs
    """
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


def fromLabramHR(filepath, **kwargs):
    data = wt.data.from_LabRAM(filepath, **kwargs)
    # add attrs
    if data.ndim == 1:
        data.attrs["dtype"] = "spectrum"
    elif data.ndim == 2:
        data.attrs["dtype"] = "spectrum" if "index" in data.variable_names else "spectralprofile"
    elif data.ndim == 3:
        data.attrs["dtype"] = "hyperspectral"

    acq_type = "Raman" if data.wl.units == "wn" else "PL"
    data.attrs['acquisition'] = f'Horiba_{acq_type}'

    return data


def fromLabramHRTimedSeries(filedir):
    raise NotImplementedError


def fromAramis(filepath):
    raise NotImplementedError

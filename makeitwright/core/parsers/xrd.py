import numpy as np
import WrightTools as wt
from ..helpers import norm

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


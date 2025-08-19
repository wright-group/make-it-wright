import numpy as np
import WrightTools as wt
from ..helpers import norm

def fromSP130(fpath, name=None):
    if fpath.split('.')[-1] != 'asc':
        print(f"filetype .{fpath.split('.')[-1]} not supported")
    else:
        with open(fpath) as f:
            txt = f.readlines()
        header_size = 0
        for i, line in enumerate(txt):
            if 'Title' in line.split() and name is None:
                name = line.split()[-1]
            if '*BLOCK' in line:
                header_size = i+1

    arr = np.genfromtxt(fpath, delimiter=',', skip_header=header_size, skip_footer=1)
    t = arr[:,0]
    sig = arr[:,1]
    t = t-t[np.argmax(sig)]

    out = wt.Data(name=name)
    out.create_variable('t', values=t, units='ns')
    out['t'].attrs['label'] = "time (ns)"
    out.create_channel('sig', values=sig)
    out['sig'].attrs['label'] = "PL counts"
    out.transform('t')
    out.create_channel('norm', values=norm(out['sig'][:], 0.01, 1))
    out['norm'].attrs['label'] = "norm. PL counts"
    
    return out


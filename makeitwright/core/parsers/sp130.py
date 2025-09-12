import WrightTools as wt


def fromSP130(fpath, **kwargs):
    data = wt.data.from_spcm(fpath, **kwargs)
    data.rename_variables(time="t")
    data.rename_channels(counts="sig")
    data.create_variable("t_shifted", values=data.t[:] - data.t[data.sig[:].argmax()], units=data.t.units)
    
    return data

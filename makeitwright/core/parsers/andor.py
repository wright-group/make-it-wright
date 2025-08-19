import WrightTools as wt
import numpy as np
import pathlib
from os import fspath


def fromAndorNeo(fpath, name=None, objective_lens='prompt', cps=False):
    """Create a data object from Andor Solis software (ascii exports).

    Parameters
    ----------
    fpath : path-like
        Path to file (should be .asc format).
        Can be either a local or remote file (http/ftp).
        Can be compressed with gz/bz2, decompression based on file name.
    name : string (optional)
        Name to give to the created data object. If None, filename is used.
        Default is None.

    Returns
    -------
    data
        New data object.
    """

    objective_lenses = {
        '5x-Jin' : 0.893,
        '20x-Jin' : 3.52, 
        '100x-Wright' : 18.2,
        '5' : 0.893,
        '20' : 3.52,
        '100' : 18.2,
        5 : 0.893,
        20 : 3.52,
        100 : 18.2
    }

    # parse filepath
    filepath = pathlib.Path(fpath)

    if not ".asc" in filepath.suffixes:
        wt.exceptions.WrongFileTypeWarning.warn(filepath, ".asc")
    # parse name
    if name is None:
        name = filepath.name.split("/")[-1]

    if objective_lens=='prompt':
        objective_lens = input(f'enter magnification for data at {name}: ')
        if not objective_lens:
            objective_lens = 0

    # create data
    ds = np.DataSource(None)
    f = ds.open(fspath(fpath), "rt")
    axis0 = []
    arr = []
    attrs = {}

    line0 = f.readline().strip()[:-1]
    line0 = [float(x) for x in line0.split(",")]  # TODO: robust to space, tab, comma
    axis0.append(line0.pop(0))
    arr.append(line0)
    
    def get_frames(f, arr, axis0):
        axis0_written = False
        while True:
            line = f.readline().strip()[:-1]
            if len(line) == 0:
                break
            else:
                line = [float(x) for x in line.split(",")]
                # signature of new frames is restart of axis0
                if not axis0_written and (line[0] == axis0[0]):
                    axis0_written = True
                if axis0_written:
                    line.pop(0)
                else:
                    axis0.append(line.pop(0))
                arr.append(line)
        return arr, axis0

    arr, axis0 = get_frames(f, arr, axis0)
    nframes = len(arr) // len(axis0)

    i = 0
    while i < 3:
        line = f.readline().strip()
        if len(line) == 0:
            i += 1
        else:
            try:
                key, val = line.split(":", 1)
            except ValueError:
                pass
            else:
                attrs[key.strip()] = val.strip()

    f.close()

    #create data object
    arr = np.array(arr)
    axis0 = np.array(axis0)
    data = wt.Data(name=name)
    if float(attrs["Grating Groove Density (l/mm)"]) == 0:
        xname = 'x'
        dtype = 'image'
        try:
            axis0 = axis0/objective_lenses[objective_lens]
            xunits = 'µm'
        except KeyError:
            xunits = 'px'
    else:
        xname = 'wl'
        xunits = 'nm'
        dtype = 'spectralprofile'
        
    axis1 = np.arange(arr.shape[-1])
    yname='y'
    try:
        axis1 = axis1/objective_lenses[objective_lens]
        yunits = 'µm'
    except KeyError:
        yunits = 'px'
    
    axes = [xname, yname]

    if nframes == 1:
        arr = np.array(arr)
        data.create_variable(name=xname, values=axis0[:, None], units=xunits)
        data.create_variable(name=yname, values=axis1[None, :], units=yunits)
    else:
        frames = np.arange(nframes)
        try:
            ct = float(attrs["Kinetic Cycle Time (secs)"])
            frames = frames*ct
            tname = 't'
            tunits = 's'
        except KeyError:
            tname = 'frame'
            tunits = None
        arr = np.array(arr).reshape(nframes, len(axis0), len(arr[0]))
        data.create_variable(name=tname, values=frames[:, None, None], units=tunits)
        data.create_variable(name=xname, values=axis0[None, :, None], units=xunits)
        data.create_variable(name=yname, values=axis1[None, None, :], units=yunits)
        axes = [tname] + axes

    if xname=='wl':
        if xunits=='nm':
            data[xname].attrs['label'] = "wavelength (nm)"
        if xunits=='wn':
            data[xname].attrs['label'] = "wavenumber (cm-1)"
    if xname=='x':
        data[xname].attrs['label'] = "x (µm)"
    if yname=='y':
        data[yname].attrs['label'] = "y (µm)"

    data.transform(*axes)
    if cps:
        try:
            arr = arr/float(attrs["Exposure Time (secs)"])
        except KeyError:
            pass
        try:
            arr = arr/int(attrs["Number of Accumulations"])
        except KeyError:
            pass

    data.create_channel(name='sig', values=arr, signed=False)
    if cps:
        data['sig'].attrs['label'] = "intensity (cps)"
    else:
        data['sig'].attrs['label'] = "counts"

    for key, val in attrs.items():
        data.attrs[key] = val

    # finish
    print("data created at {0}".format(data.fullpath))
    print("  axes: {0}".format(data.axis_names))
    print("  shape: {0}".format(data.shape))
    data.attrs['dtype']=dtype

    return data


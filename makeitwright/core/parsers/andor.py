import WrightTools as wt


# builtin conversion from pixels to um (from somebody's records)
# for different objectives
# roughly, pixel size in microns / magnification
px_per_um = {
    '5x-Jin' : 0.893,
    '20x-Jin' : 3.52, 
    '100x-Wright' : 18.2,
    '5' : 0.893,
    '20' : 3.52,
    '100' : 18.2,
}


def fromAndorNeo(fpath, name=None, px_per_um=None):
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
    px_per_um : float-like (optional)
        if present, camera spatial dimensions will be mapped in micron units.
        if not present, spatial variables of camera will be a unitless index

    Returns
    -------
    data
        New data object.
    """
    # parse filepath
    data:wt.Data = wt.data.from_Solis(fpath, name=name, verbose=True)
    data.rename_variables(xindex="x", yindex="y", wm="wl")

    for var in (("x", "y") and data.variable_names):
        if px_per_um:
            data[var] /= px_per_um
            data[var].units = 'µm'

    dtype = "image" if "x" in data.variable_names else "spectralprofile"
    data.attrs.update(dtype=dtype)

    if "wl" in data.variable_names:
        data["wl"].attrs['label'] = "wavelength (nm)" if data["wl"].units == "nm" else "wavenumber (cm-1)"

    if data.signal.units == "Hz"
        data.signal.label = "intensity (cps)"
    else:
        data.signal.label = "counts"

    return data


import pathlib
import os
import warnings
import numpy as np
import WrightTools as wt

from . import image
from . import spectralprofile
from .helpers import roi, set_label, get_channels
import makeitwright.styles as styles

APD_PIXEL = (1325, 1080)
SLIT_PIXEL_COLUMN = 1325
CMOS_RESPONSE_CUTOFF_NM = 370
MATERIAL_OPTICAL_DATA = { #R or 1-T vs nm
    'UVFS':'transmittance_references/UVFS.csv', #UV fused silica
    'sapphire':'transmittance_references/sapphire.csv',
    'soda lime':'transmittance_references/soda-lime.csv', #soda lime glass
    'ThorLabs silver':'transmittance_references/Ag-P01.csv', #ThorLabs -P01 silver
    'BK7':'transmittance_references/BK7.csv', #Schott BK7 borosilicate glass
    'fluorite':'transmittance_references/CaF2.csv',
    'MgF2':'transmittance_references/MgF2.csv',
    'ThorLabs UV Al':'transmittance_references/UVAl.csv' #ThorLabs -F01 UV-enhanced aluminum
    }

def plot_image(data, channel=0, **kwargs):  
    params = {}
    params.update(styles.image_andor)
    params.update(**kwargs)
    image.plot_image(data, channel, **params)

def plot_profile(data, channel=0, **kwargs):
    params = {}
    params.update(styles.profile_andor)
    params.update(**kwargs)  
    spectralprofile.plot_profile(data, channel, **params)
    
def plot_decomposition(data, channel=0, mode='R', **kwargs):
    warnings.warn("plot_decomposition will be removed from the next version of makeitwright. Its functionality will be merged into the method plot1D of the artists module.", DeprecationWarning)

    modes = {
        'A' : styles.decomposition_andor_A,
        'PL' : styles.decomposition_andor_PL,
        'R' : styles.decomposition_andor_R,
        'RR0' : styles.decomposition_andor_RR0,
        'T' : styles.decomposition_andor_T
        }
     
    params = {}
    params.update(modes[mode])
    params.update(**kwargs)
   
    spectralprofile.plot_decomposition(data, 'wl', 'y', channel, **params)
        
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
    filestr = os.fspath(fpath)
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
    f = ds.open(filestr, "rt")
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

def _get_reference_material_array(data, channel, material):
    """
    Return a transmittance array of a material interpolated from data wavelength values.
    """
    if material not in MATERIAL_OPTICAL_DATA:
        print(f"Material {material} was not found in the reference data directory.")
        return 1
    ch = get_channels(data, channel)[0]
    data_wavelengths = data.axes[0].points
    material_array = np.genfromtxt(MATERIAL_OPTICAL_DATA[material], delimiter=',', skip_header=1)
    material_wavelengths = material_array[:,0].flatten()
    material_transmittance = material_array[:,1].flatten()
    reference_array = np.interp(data_wavelengths, material_wavelengths, material_transmittance).reshape((data_wavelengths.size, 1))
    reference_array = np.repeat(reference_array, data[ch].shape[1], axis=1)
    return reference_array

def _get_processed_signal_arrays(sample_data, reference_data, 
                            sample_channel, reference_channel,
                            dark_reference_data, dark_reference_channel,
                            dark_wavelength_range,
                            background_ROI):
    #extract channels
    ch_sample, ch_reference = get_channels(sample_data, sample_channel)[0], get_channels(reference_data, reference_channel)[0]
    signal_sample, signal_reference = sample_data[ch_sample][:], reference_data[ch_reference][:]
    #apply dark count correction
    if dark_reference_data is not None:
        ch_dark = get_channels(dark_reference_data, dark_reference_channel)[0]
        dark = dark_reference_data[ch_dark][:]
        signal_sample -= dark
        signal_reference -= dark
    elif dark_wavelength_range is not None:
        dark = get_spectral_dark_counts(sample_data, sample_channel, dark_wavelength_range)
        signal_sample -= dark
        dark = get_spectral_dark_counts(reference_data, reference_channel, dark_wavelength_range)
        signal_reference -= dark
    #convert data to cps
    if sample_data[ch_sample].units != 'Hz':
        signal_sample /= get_exposure_time(sample_data)
    if reference_data[ch_reference].units != 'Hz':
        signal_reference /= get_exposure_time(reference_data)
    #identify correction coefficient from y background
    if background_ROI is not None:
        background_sample = roi(sample_data, background_ROI, return_arrs=True)[ch_sample]
        background_reference = roi(reference_data, background_ROI, return_arrs=True)[ch_reference]
        correction_factor = np.average(background_sample/background_reference)
        signal_reference *= correction_factor
    signals = {
        'sample' : signal_sample,
        'reference' : signal_reference
    }
    return signals

def get_APD_location(data):
    """
    Return the (x,y) coordinates corresponding to the typical location of the APD single photon counter.
    """
    return image.get_pixel_location(data, APD_PIXEL)

def get_slit_location(data):
    """
    Return the x-axis coordinate corresponding to the nominal center of the spectrograph entrance slit.
    """
    return data.axes[0].points[SLIT_PIXEL_COLUMN]

def get_exposure_time(data):
    """
    Return the total exposure time of an Andor Neo acquisition.
    """
    try:
        exposure = float(data.attrs["Exposure Time (secs)"])
    except KeyError:
        exposure = 1
    try:
        accumulations = int(data.attrs["Number of Accumulations"])
    except KeyError:
        accumulations = 1
    total_time = exposure*accumulations
    return total_time

def get_background_spectrum(data, channel, background_y_range):
    """
    Return a background spectrum array by averaging along a background y-range in a hyperspectral profile.
    """
    ch = get_channels(data, channel)[0]
    background = roi(data, {1:(list(background_y_range), 'average')}, return_arrs=True)[ch]
    return background

def get_background_y_profile(data, channel, background_wavelength_range):
    """
    Return a background y-profile by averaging along a background wavelength range in a hyperspectral profile.
    """
    ch = get_channels(data, channel)[0]
    background = roi(data, {0:(list(background_wavelength_range), 'average')}, return_arrs=True)[ch]
    return background

def get_image_dark_counts(data, channel, slit_edges):
    """
    Return average dark counts using regions of an image obscured by slit blades as a reference.
    """
    ch = get_channels(data, channel)[0]
    ROI_left = roi(data, {0:[0,slit_edges[0]]}, return_arrs=True)
    ROI_right = roi(data, {0:[slit_edges[1]]}, return_arrs=True)
    dark_counts = (np.average(ROI_left[ch])+np.average(ROI_right[ch]))/2
    return dark_counts

def get_spectral_dark_counts(data, channel, dark_wavelength_range):
    """
    Return average dark counts using a wavelength range approximated to be dark as a reference.
    """
    if data.axes[0].min()>dark_wavelength_range[1]:
        raise ValueError(f"The minimum data wavelength {data.axes[0].min()} is beyond the maximum dark wavelength value {dark_wavelength_range[1]}")
    ch = get_channels(data, channel)[0]
    ROI = roi(data, {0:list(dark_wavelength_range)}, return_arrs=True)
    dark_counts = np.average(ROI[ch])
    return dark_counts

def compute_cps(data, channel, name=None):
    """
    Generate a new channel from an existing one in units of counts per second.
    """
    ch = get_channels(data, channel)[0]
    if data[ch].units == 'Hz':
        print(f"Channel {ch} already has units of cps.")
    else:
        if name is None:
            name = f'{ch}_cps'
        exposure_time = get_exposure_time(data)
        data.create_channel(name, values=data[ch][:]/exposure_time, units='Hz')
        if data[ch].signed:
            data[name].signed=True        

def compute_reflectance(sample_data, reference_data, 
                            sample_channel=0, reference_channel=0,
                            dark_reference_data=None, dark_reference_channel=0,
                            dark_wavelength_range=[0, CMOS_RESPONSE_CUTOFF_NM],
                            background_ROI=None,
                            reference_material='BK7'):
    """
    Determine the spectral reflectance of a sample using the spectral reflectance of the substrate underneath.
    Optional arguments consist of various corrections to improve the quantitative accuracy of the spectrum.
    This method is designed for hyperspectral y-traces with dimensions (wavelength, y), corresponding to a spectrally dispersed y-axis of the region imaged by the camera.

    Arguments
    ---------
    sample_data : WrightTools.Data - The spectral reflection of the sample.
    reference_data : WrightTools.Data - The spectral reflection of the substrate.
    sample_channel, reference_channel : str or int, optional - The keys or indices of the data channels to compute reflectance from. Default is 0, the raw intensities.
    dark_reference_data : WrightTools.Data or None, optional - A camera frame taken without illumination for dark count correction. Default is None.
        Note that the dark reference should be collected under the same exposure conditions as the sample and reference data, otherwise dark count correction will be inaccurate.
    dark_reference_channel : str or int, optional - The channel corresponding to the dark frame. Default is 0, the raw intensity.
    dark_wavelength_range : iterable of 2 numbers or None, optional - The range of wavelengths that can be referenced for an approximation of dark counts and applied for dark count correction. 
        Default range is from 0 to a pre-determined responsivity cutoff of the instrument.
    background_ROI : dict or None, optional - A region of the data that contains only substrate in both the sample and reference. Dictionary is formatted according to processhelpers.roi()
        Default is None. If a region is provided, background correction will be implemented under the assumption that the reflectivity of the substrate should be identical in both frames.     
    reference_material : str or None, optional - Key for substrate reference data to convert relative reflectance into absolute reflectance.
        Default is BK7 for N-BK7 borosilicate optical glass. For available data and corresponding keys, see andor.MATERIAL_OPTICAL_DATA.
        If None, the reflectance is computed relative to the substrate and labelled as such.
    
    Returns
    -------
    None - Creates a new channel in sample_data that contains the computed reflectance frame.
    """
    #apply corrections to signal arrays
    signals = _get_processed_signal_arrays(sample_data, reference_data, sample_channel, reference_channel, dark_reference_data, dark_reference_channel, dark_wavelength_range, background_ROI)
    r, r0 = signals['sample'], signals['reference']
    #compute the reflectance array
    R = r/r0
    #mask obviously impossible values
    R = np.where(R>100, 100, R)
    R = np.where(R<-100, -100, R)
    #correct to absolute reflectance if possible
    name = 'R_rel'
    label = "reflectance vs. substrate"
    if reference_material is not None:
        substrate_reflectance = _get_reference_material_array(sample_data, sample_channel, reference_material)
        if isinstance(substrate_reflectance, np.ndarray):
            R *= substrate_reflectance
            name = 'R'
            label = "reflectance"
    #add channel to sample data
    sample_data.create_channel(name, values=R)
    set_label(sample_data, name, label)

def compute_reflection_contrast(sample_data, reference_data,
                            sample_channel=0, reference_channel=0,
                            dark_reference_data=None, dark_reference_channel=0,
                            dark_wavelength_range=[0, CMOS_RESPONSE_CUTOFF_NM],
                            background_ROI=None):
    """
    Determine the spectral reflection contrast of a sample vs. the substrate underneath.
    Optional arguments consist of various corrections to improve the quantitative accuracy of the spectrum.
    This method is designed for hyperspectral y-traces with dimensions (wavelength, y), corresponding to a spectrally dispersed y-axis of the region imaged by the camera.

    Arguments
    ---------
    sample_data : WrightTools.Data - The spectral reflection of the sample.
    reference_data : WrightTools.Data - The spectral reflection of the substrate.
    sample_channel, reference_channel : str or int, optional - The keys or indices of the data channels to compute reflection contrast from. Default is 0, the raw intensities.
    dark_reference_data : WrightTools.Data or None, optional - A camera frame taken without illumination for dark count correction. Default is None.
        Note that the dark reference should be collected under the same exposure conditions as the sample and reference data, otherwise dark count correction will be inaccurate.
    dark_reference_channel : str or int, optional - The channel corresponding to the dark frame. Default is 0, the raw intensity.
    dark_wavelength_range : iterable of 2 numbers or None, optional - The range of wavelengths that can be referenced for an approximation of dark counts and applied for dark count correction. 
        Default range is from 0 to a pre-determined responsivity cutoff of the instrument.
    background_ROI : dict or None, optional - A region of the data that contains only substrate in both the sample and reference. Dictionary is formatted according to processhelpers.roi()
        Default is None. If a region is provided, background correction will be implemented under the assumption that the reflectivity of the substrate should be identical in both frames.     

    Returns
    -------
    None - Creates a new channel in sample_data that contains the computed reflection contrast frame.
    """
    #apply corrections to signal arrays
    signals = _get_processed_signal_arrays(sample_data, reference_data, sample_channel, reference_channel, dark_reference_data, dark_reference_channel, dark_wavelength_range, background_ROI)
    r, r0 = signals['sample'], signals['reference']
    #compute the reflectance contrast array
    RR0 = (r-r0)/r0
    #mask obviously impossible values
    RR0 = np.where(RR0>100, 100, RR0)
    RR0 = np.where(RR0<-100, -100, RR0)
    #add channel to sample data
    sample_data.create_channel('RR0', values=RR0)
    sample_data['RR0'].signed = True
    set_label(sample_data, 'RR0', "reflection contrast")

def compute_transmittance(sample_data, reference_data, 
                            sample_channel=0, reference_channel=0,
                            dark_reference_data=None, dark_reference_channel=0,
                            dark_wavelength_range=[0, CMOS_RESPONSE_CUTOFF_NM],
                            background_ROI=None):
    """
    Determine the spectral transmittance of a sample using the spectral transmittance of the substrate underneath.
    Optional arguments consist of various corrections to improve the quantitative accuracy of the spectrum.
    This method is designed for hyperspectral y-traces with dimensions (wavelength, y), corresponding to a spectrally dispersed y-axis of the region imaged by the camera.

    Arguments
    ---------
    sample_data : WrightTools.Data - The spectral transmission of the sample.
    reference_data : WrightTools.Data - The spectral transmission of the substrate.
    sample_channel, reference_channel : str or int, optional - The keys or indices of the data channels to compute transmittance from. Default is 0, the raw intensities.
    dark_reference_data : WrightTools.Data or None, optional - A camera frame taken without illumination for dark count correction. Default is None.
        Note that the dark reference should be collected under the same exposure conditions as the sample and reference data, otherwise dark count correction will be inaccurate.
    dark_reference_channel : str or int, optional - The channel corresponding to the dark frame. Default is 0, the raw intensity.
    dark_wavelength_range : iterable of 2 numbers or None, optional - The range of wavelengths that can be referenced for an approximation of dark counts and applied for dark count correction. 
        Default range is from 0 to a pre-determined responsivity cutoff of the instrument.
    background_ROI : dict or None, optional - A region of the data that contains only substrate in both the sample and reference. Dictionary is formatted according to processhelpers.roi()
        Default is None. If a region is provided, background correction will be implemented under the assumption that the transmissivity of the substrate should be identical in both frames.     
    reference_material : str or None, optional - Key for substrate reference data to convert relative transmittance into absolute transmittance.

    Returns
    -------
    None - Creates a new channel in sample_data that contains the computed transmittance frame.
    """
    #apply corrections to signal arrays
    signals = _get_processed_signal_arrays(sample_data, reference_data, sample_channel, reference_channel, dark_reference_data, dark_reference_channel, dark_wavelength_range, background_ROI)
    t, t0 = signals['sample'], signals['reference']
    #compute the reflectance array
    T = t/t0
    #mask obviously impossible values
    T = np.where(T>100, 100, T)
    T = np.where(T<-100, -100, T)
    #add channel to sample data
    sample_data.create_channel('T', values=T)
    set_label(sample_data, 'T', "transmittance")

def compute_absorbance(reflectance_data, transmittance_data,
                                    reflectance_channel='R', transmittance_channel='T',
                                    name='A'):
    """
    Compute the directional absorbance of a material from its reflectance and transmittance data.
    
    Arguments
    ---------
    reflectance_data, transmittance_data : WrightTools.Data - The reflectance and transmittance data.
    reflectance_channel, transmittance_channel : str or int, optional - The keys or indices of the channels containing the spectral reflectance and transmittance. Defaults are 'R' and 'T'.
    name : str, optional - The name for the new absorbance channel. Default is 'A'.

    Returns
    -------
    None - Creates a new channel in each data object for the computed directional absorbance.
    """
    chR = get_channels(reflectance_data, reflectance_channel)[0]
    chT = get_channels(transmittance_data, transmittance_channel)[0]
    R = reflectance_data[chR][:]
    T = transmittance_data[chT][:]
    transmittance_data.create_channel(name, values=1-T-R)
    set_label(transmittance_data, name, "absorbance")
    reflectance_data.create_channel(name, values=1-T-R)
    set_label(reflectance_data, name, "absorbance")

def correct_PL_background(data, ybkg=None, cps=True):
    """
    Apply dark count corrections to PL data.
    """
    if type(data) is not list:
        data = [data]

    for d in data:
        k = 0
        if ybkg is not None:
            k = roi(d, {'y':(ybkg,'average')}, return_arrs=True)['sig']
            k = k.reshape((k.size,1))
            k = k.repeat(d['sig'][:].shape[1], axis=1)

        cycle_time = 1
        if cps:
            try:
                exposure = float(d.attrs["Exposure Time (secs)"])
            except KeyError:
                exposure = 1
            try:
                accum = int(d.attrs["Number of Accumulations"])
            except KeyError:
                accum = 1
            cycle_time = exposure*accum

        sig = d['sig'][:]
        if cps:
            d.create_channel('corr', values=(sig-k)/cycle_time, units='Hz')
            label = "PL intensity (cps)"
        else:
            d.create_channel('corr', values=(sig-k)/cycle_time)
            label = "PL counts"
        set_label(d, 'corr', label)

        if len(data)==1:
            print(f"PL correction applied to data {data[0].natural_name}")
        else:
            print(f"PL correction applied to {len(data)} spectra")
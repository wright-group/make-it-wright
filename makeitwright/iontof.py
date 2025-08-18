import numpy as np
import cmocean
from .lib import hyperspectral, styles, helpers


def relative_proportion(data, channel0, channel1):
    """
    Calculate the relative proportion of signal between two channels on a scale of -1 to 1, with 0 being equal in proportion.
    
    Parameters
    ----------
    data : Data object of WrightTools data module
        The data to which the channels of interest belong.
        
    channel0, channel1 : str or int
        The natural names or indices of the channels to compare.

    Returns
    -------
    None.
    """
    
    #convert channel indices to natural names
    channels = helpers.parse_args(data, channel0, channel1, dtype='Channel')
    
    #ensure regions with no signal don't show up
    hyperspectral.remove_background(data, *channels, threshold_value=0.99, new_value=1)
    dump = [data.channels[-2].natural_name, data.channels[-1].natural_name]
    
    #calculate normalized difference between channels
    helpers.normalize_by_axis(data, channels[0], 'x', 'y', 'scan')
    helpers.normalize_by_axis(data, channels[1], 'x', 'y', 'scan')
    ch_arr0 = data.channels[-2][:]
    dump.append(data.channels[-2].natural_name)
    ch_arr1 = data.channels[-1][:]
    dump.append(data.channels[-1].natural_name)        
    ch_arr = (ch_arr0-ch_arr1)/(ch_arr0+ch_arr1)
    ch_arr = np.nan_to_num(ch_arr)
    
    for trash in dump:
        data.remove_channel(trash, verbose=False)
    
    #create new channel
    ch_name = "rel_mass_"+channels[0]+"_"+channels[1]
    data.create_channel(ch_name, values=ch_arr, verbose=True)
    data[ch_name].signed = True

def plot_image(data, channel, **kwargs):
    
    params = {}
    params.update(styles.image_iontof)
    if data[channel].signed:
        params["cmap"] = cmocean.cm.curl
    params.update(**kwargs)

    hyperspectral.plot_image(data, channel, **params)

def plot_profile(data, profile_axis, channel, **kwargs):
    
    params = {}
    params.update(styles.profile_iontof)
    if data[channel].signed:
        kwargs["cmap"] = cmocean.cm.curl
    params.update(**kwargs)
    
    hyperspectral.plot_profile(data, profile_axis, channel, **params)

def plot_depth_trace(data, channel, **kwargs):
    
    params = {}
    params.update(styles.decomposition_iontof)
    params.update(**kwargs)
    
    hyperspectral.plot_decomposition(data, 'x', 'y', 'scan', channel, **kwargs)


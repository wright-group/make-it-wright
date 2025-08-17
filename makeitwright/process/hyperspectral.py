__name__ = "hyperspectral"
__author__ = "Chris Roy, Song Jin Research Group, Dept. of Chemistry, University of Wisconsin - Madison"

"""
Processing and plotting methods for 3-dimensional WrightTools data objects containing two spatial axes and one non-spatial axis (i.e. a spectral axis).

The two spatial axes are generally defined as x and y. 
The third axis is referred to as "spectral axis", 
but an arbitrary non-spatial or pseudo-spatial axis may be used where relevant.

Data axes must be ordered (spatial x, spatial y, non-spatial).
"""

#import
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
import processhelpers as proc
import styles

def remove_background(data, channel, threshold=0.1, negative=False, return_mask=False, max_ref_count=10):
    """
    Remove background pixels from the x-y plane of the hyperspectral image using the spectrally binned image as a reference signal.
    Background x-y points will be set to 0 along the entire spectral axis.
    
    Parameters
    ----------
    data : Data object of WrightTools data module.
        The data for which the background-subtracted channel will be generated.
    channel : int or str
        The channel that will be duplicated with subtracted background.
    threshold : float between 0 and 1, optional
        The fraction of the maximum reference value below which is to be considered background. 
        The default is 0.1.
    negative : bool, optional
        Subtract everything but the background instead. Useful if the region of interest is the signal minimum. 
        The default is False.
    max_ref_count : int, optional
        The number of highest x-y points in the image to be averaged as a reference for the channel maximum. Useful to avoid false maxima caused by spikes.
        The default is 10.

    Returns
    -------
    None.
        Adds a background-subtracted channel to the Data instance.
    """
    #parse channel argument as str
    channel, = proc.parse_args(data, channel, dtype='Channel')
    
    #generate a spectrally binned image of the channel
    ch_arr = np.sum(data[channel][:], axis=2)
    #get a representative maximum value from the reference image
    if max_ref_count > ch_arr.size:
        max_ref_count = int(ch_arr.size/10)
    ordered = np.sort(ch_arr.flatten())[-max_ref_count:]
    ch_max = np.average(ordered)
    #generate a mask array
    bkg = np.where(ch_arr < threshold*ch_max, 0, ch_arr)
    bkg = np.where(bkg > 0, 1, bkg)
    #remove background using mask
    mask = np.repeat(bkg[:,:,None], data.axes[2].size, axis=2)
    if negative:
        mask = 1-mask
    nobkg = data[channel][:] * mask
    #create background-subtracted channel
    data.create_channel(name=channel+"_nobkg", values=nobkg)
    data[channel+'_nobkg'].signed = data[channel].signed
    if return_mask:
        data.create_channel(name=data[channel].natural_name+"_mask", values=mask)
    
def get_profile(data, profile_axis, ROI=None):
    """
    Extract profile from an arbitrary pair of points in a selected 2D subspace of the data's 3 dimensions.
    Arguments
    ------------------------------------
    data: WrightTools Data instance. Must have at least 3 axes.
    ROI should be a dict object containing the beginning and end points for each relevant axis, returns data object
    """
    profile_axis, = proc.parse_args(data, profile_axis)
    non_profile_axis = [axis.natural_name for axis in data.axes if axis.natural_name != profile_axis][0]
    spectral_axis = [axis.natural_name for axis in data.axes if axis.natural_name != profile_axis][1]

    dims_too_low = False
    if ROI is not None:
        if [val for val in ROI.values() if val == 'all']:
            dims_too_low = True
        if spectral_axis in ROI.keys():
            if type(ROI[spectral_axis]) is int or type(ROI[spectral_axis]) is float:
                dims_too_low = True
        if profile_axis in ROI.keys():
            if type(ROI[profile_axis]) is int or type(ROI[profile_axis]) is float:
                dims_too_low = True
        if len([val for val in ROI.values() if type(val) is int or type(val) is float])>1:
            dims_too_low = True
        if dims_too_low:
            print("Dimensionality of ROI is too low. Do not collapse any dimensions of the data before calling this method.")
            return
        if not dims_too_low:
            out = proc.roi(data, ROI)
    else:
        out = data
    if len(out.axes) > 2:
        out = proc.roi(out, {non_profile_axis:'all'})
    
    out.transform()
    
    for channel in out.channels:
        ch_name = channel.natural_name
        ch_values = out[ch_name][:].transpose()
        out.remove_channel(ch_name, verbose=False)
        out.create_channel(ch_name, values=ch_values, verbose=False)
    for variable in out.variables:
        var_name = variable.natural_name
        var_values = out[var_name][:].transpose()
        var_units = variable.units
        out.remove_variable(var_name, verbose=False)
        out.create_variable(var_name, values=var_values, units=var_units, verbose=False)
    
    out.transform(spectral_axis, profile_axis)
    print(f'profile along direction <{profile_axis}> extracted')
    
    return out
    
def plot_image(data, channel, **kwargs):
    
    #convert axis/channel indices to natural names
    channel, = proc.parse_args(data, channel, dtype='Channel')
    non_spatial_axis = data.axes[-1].natural_name

    #set parameters for plotting from kwargs
    params = {                            
        "ROI" : None,
        'ticks' : 'auto',
        "vrange" : None,
        "title" : None
        }
    params.update(styles.image)
    if data[channel].signed:
        params["cmap"] = mpl.cm.RdBu_r
    params.update(**kwargs)
    
    #extract ROI
    if params["ROI"] is not None:
        out = proc.roi(data, params["ROI"])
        if len(out.axes) != 2:
            out = proc.roi(out, {non_spatial_axis : 'all'})
    else:
        out = proc.roi(data, {non_spatial_axis : 'all'})
    
    #determine range to be plotted
    if params["vrange"] is None:
        vrange = proc.get_range(out, reference_key=channel)
    else:
        vrange = params["vrange"]
    
    #setup plot frame
    fig, ax = plt.subplots(figsize=(params['fig_width'], params['fig_height']))
    
    #plot data
    xgrid, ygrid = np.meshgrid(out.axes[0][:], out.axes[1][:])
    mesh = ax.pcolormesh(xgrid, ygrid, np.transpose(out[channel][:]), cmap=params["cmap"], vmin=vrange[0], vmax=vrange[1])
    ax.set_xlabel(params["xlabel"])
    ax.set_ylabel(params["ylabel"])
    ax.set_aspect("equal")
    if params["title"] is not None:
        ax.set_title(params["title"])
    
    #set ticks
    if params['ticks'] == 'auto':
        ticks = np.linspace(vrange[0], vrange[1], num=11)
    elif params['ticks'] is None:
        ticks = []
    else:
        ticks = params['ticks']
    
    # plot colorbar
    cbar = plt.colorbar(mesh)
    cbar.set_ticks(ticks)
    cbar.set_label(params["cbar_label"])

def plot_profile(data, profile_axis, channel, **kwargs):
    
    #convert axis/channel indices to natural names
    profile_axis, = proc.parse_args(data, profile_axis)
    spectral_axis = data.axes[-1].natural_name
    channel, = proc.parse_args(data, channel, dtype='Channel')
    non_profile_axis = [axis.natural_name for axis in data.axes[:-1] if axis.natural_name != profile_axis][0]

    #set parameters for plotting from kwargs
    params = {
        "ROI" : None,
        'ticks' : 'auto',
        "vrange" : None,
        "reference_lines" : None,
        "title" : None
        }
    params.update(styles.profile)
    if data[channel].signed:
        params["cmap"] = mpl.cm.RdBu_r
    params.update(**kwargs)
    
    #extract ROI
    if params["ROI"] is not None:
        out = proc.roi(data, params["ROI"])
        if len(out.axes) != 2:
            out = proc.roi(out, {non_profile_axis : 'all'})
    else:
        out = proc.roi(data, {non_profile_axis : 'all'})
    out.transform(spectral_axis, profile_axis)

    #determine range to be plotted
    if params["vrange"] is None:
        vrange = proc.get_range(out, reference_key=channel)
    else:
        vrange = params["vrange"]
    
    #setup plot frame
    fig, ax = plt.subplots(figsize=(params['fig_width'], params['fig_height']))
    
    #plot data
    xgrid, ygrid = np.meshgrid(out.axes[1][:], out.axes[0][:])
    try:
        mesh = ax.pcolormesh(xgrid, ygrid, np.transpose(out[channel][:]), cmap=params["cmap"], vmin=vrange[0], vmax=vrange[1])
    except TypeError:
        mesh = ax.pcolormesh(xgrid, ygrid, out[channel][:], cmap=params["cmap"], vmin=vrange[0], vmax=vrange[1])
    ax.set_xlabel(params["xlabel"])
    ax.set_ylabel(params["ylabel"])
    
    if params["reference_lines"] is not None:
        if type(params["reference_lines"]) is not list:
            params["reference_lines"] = [params["reference_lines"]]
        for reference_line in params["reference_lines"]:
            ax.axvline(x=reference_line, linewidth=1, color='grey', linestyle='--', alpha=0.25)
    
    #set ticks
    if params['ticks'] == 'auto':
        ticks = np.linspace(vrange[0], vrange[1], num=11)
    elif params['ticks'] is None:
        ticks = []
    else:
        ticks = params['ticks']
    
    # plot colorbar
    cbar = plt.colorbar(mesh)
    cbar.set_ticks(ticks)
    cbar.set_label(params["cbar_label"])
    
    if params["title"] is not None:
        ax.set_title(params["title"])
    
def plot_decomposition(data, x_axis, y_axis, spectral_axis, channel, **kwargs):
    #convert axis/channel indices to natural names
    x_axis, y_axis, spectral_axis = proc.parse_args(data, x_axis, y_axis, spectral_axis)
    channel, = proc.parse_args(data, channel, dtype='Channel')

    #set parameters for plotting from kwargs
    params = {
        "ROI" : None,
        "xrange" : None,
        "vrange" : None,
        "yscale" : 'linear',
        "binning" : None,
        "reference_lines" : None,
        "xticks" : True,
        "yticks" : True,
        "title" : None
        }
    params.update(styles.decomposition)
    params.update(**kwargs)

    #extract ROI
    if params["ROI"] is not None:
        out = proc.roi(data, params["ROI"])
    else:
        out = data
      
    #identify spatial ranges for indexing
    xrange = proc.get_range(out, reference_key=spectral_axis, dtype='Axis')

    #setup plot frame
    fig, ax = plt.subplots(figsize=(params['fig_width'], params['fig_height']))
    
    #plot data
    if params["binning"] is not None:
        if params["binning"] == 'average':
            arr_out = np.sum(out[channel][:], axis=(0,1))/(np.count_nonzero(out[channel][:])/out[channel].shape[2])
        if params["binning"] == 'sum':
            arr_out = np.sum(out[channel][:], axis=(0,1))
        #determine range to be plotted
        vrange = proc.vrange(arr_out, out[channel].signed, window=1)

        ax.plot(out[spectral_axis].points, arr_out,
            params["marker"], linewidth=params["linewidth"], alpha=1, color=params["color"])
    else:
        #determine range to be plotted
        vrange = proc.get_range(out, reference_key=channel)
        for i in range(out[x_axis].size):
            for j in range(out[y_axis].size):
                if np.sum(out[channel][i,j,:]) != 0:
                    ax.plot(out[spectral_axis].points, out[channel][i,j,:], 
                        params["marker"], linewidth=params["linewidth"], alpha=params["alpha"], color=params["color"])
    
    if params["xrange"] is not None:
        xrange = params["xrange"]
    ax.set_xlim(*xrange)
    if not params["xticks"]:
        ax.set_xticks([])
    if params["vrange"] is not None:
        vrange = params["vrange"]
    ax.set_ylim(*vrange)
    ax.set_yscale(params["yscale"])
    if not params["yticks"]:
        ax.set_yticks([])
    
    if out[channel].signed:
        ax.axhline(y=0, color='black', linewidth=1)
    
    if params["reference_lines"] is not None:
        if type(params["reference_lines"]) is not list:
            params["reference_lines"] = [params["reference_lines"]]
        for reference_line in params["reference_lines"]:
            ax.axvline(x=reference_line, zorder=0, linewidth=1, color='grey', linestyle='--', alpha=0.5)
    
    # label plot
    ax.set_xlabel(params["xlabel"])
    ax.set_ylabel(params["ylabel"])
    ax.set_yscale(params["yscale"])
    if params["title"] is not None:
        ax.set_title(params["title"])
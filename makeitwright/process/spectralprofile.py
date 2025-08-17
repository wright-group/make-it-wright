"""Methods for two-dimensional data consisting of a spectral axis (0) and a spatial axis (1)."""


import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from . import helpers as ph
import makeitwright.styles as styles


def remove_spectral_background(data, channel, spatial_reference_range, name=None, create_background_channel=False, talkback=True):
    """
    Remove background along the spatial axis using a specified range along the other axis as reference.
    Creates a new channel with the background-subtracted array.

    Arguments
    ---------
    data : WrightTools.Data - The data.
    background_axis : str or int - The axis along which the background is.
    channel : str or int - The channel to subtract the background from.
    reference_axis_range : list or int
    
    Returns
    -------
    None - Creates new background-subtracted Channels in the Data instance. 
    """
    #identify channel and categorize axes
    channel = ph.get_channels(data, channel)[0]
    spectral_axis = ph.get_axes(data, 0)[0]
    spatial_axis = ph.get_axes(data, 1)[0]
    #construct the background array
    if isinstance(spatial_reference_range, int):
        spectral_background = ph.roi(data, {spatial_axis:spatial_reference_range}, return_arrs=True)[channel].reshape(data[spectral_axis].shape)
    else:
        spectral_background = ph.roi(data, {spatial_axis:(spatial_reference_range,'average')}, return_arrs=True)[channel].reshape(data[spectral_axis].shape)
    spatial_points = np.ones(data[spatial_axis].shape)
    background = spectral_background*spatial_points
    #create background-subtracted channel
    if name is None:
        name = f"{channel}_bkgsub_spectral"
    data.create_channel(name, values=data[channel][:]-background, units=data[channel].units)
    if data[channel].signed:
        data[name].signed = True
    if create_background_channel:
        data.create_channel(f"spectral_bkg_{channel}", values=background, units=data[channel.units])
        if data[channel].signed:
            data[name].signed = True
    
    if talkback:
        print(f"subtracted spectral background from data {data.natural_name}")

def remove_spatial_background(data, channel, spectral_reference_range, name=None, create_background_channel=False, talkback=True):
    """
    Remove background along the spatial axis using a specified range along the other axis as reference.
    Creates a new channel with the background-subtracted array.

    Arguments
    ---------
    data : WrightTools.Data - The data.
    background_axis : str or int - The axis along which the background is.
    channel : str or int - The channel to subtract the background from.
    reference_axis_range : list or int
    
    Returns
    -------
    None - Creates new background-subtracted Channels in the Data instance. 
    """
    #identify channel and categorize axes
    channel = ph.get_channels(data, channel)[0]
    spectral_axis = ph.get_axes(data, 0)[0]
    spatial_axis = ph.get_axes(data, 1)[0]
    #construct the background array
    if isinstance(spectral_reference_range, int):
        spatial_background = ph.roi(data, {spectral_axis:spectral_reference_range}, return_arrs=True)[channel].reshape(data[spatial_axis].shape)
    else:
        spatial_background = ph.roi(data, {spectral_axis:(spectral_reference_range,'average')}, return_arrs=True)[channel].reshape(data[spatial_axis].shape)
    spectral_points = np.ones(data[spectral_axis].shape)
    background = spatial_background*spectral_points
    #create background-subtracted channel
    if name is None:
        name = f"{channel}_bkgsub_spatial"
    data.create_channel(name, values=data[channel][:]-background, units=data[channel].units)
    if data[channel].signed:
        data[name].signed = True
    if create_background_channel:
        data.create_channel(f"spatial_bkg_{channel}", values=background, units=data[channel.units])
        if data[channel].signed:
            data[name].signed = True
    
    if talkback:
        print(f"subtracted spatial background from data {data.natural_name}")

def remove_combined_background(data, channel, spectral_reference_range, spatial_reference_range, name=None, create_background_channel=False, talkback=True):
    """
    Remove background from data using a range of the spectral profile along each axis as reference. The background is a matrix product of the two background arrays.
    """
    def __at(arr, val):
        return (np.abs(arr-val)).argmin()
    #identify channel and categorize axes
    channel = ph.get_channels(data, channel)[0]
    spectral_axis = ph.get_axes(data, 0)[0]
    spatial_axis = ph.get_axes(data, 1)[0]
    #extract background along each axis
    if isinstance(spatial_reference_range, int):
        spectral_background = ph.roi(data, {spatial_axis:spatial_reference_range}, return_arrs=True)[channel].reshape(data[spectral_axis].shape)
    else:
        spectral_background = ph.roi(data, {spatial_axis:(spatial_reference_range,'average')}, return_arrs=True)[channel].reshape(data[spectral_axis].shape)
    if isinstance(spectral_reference_range, int):
        spatial_background = ph.roi(data, {spectral_axis:spectral_reference_range}, return_arrs=True)[channel].reshape(data[spatial_axis].shape)
    else:
        spatial_background = ph.roi(data, {spectral_axis:(spectral_reference_range,'average')}, return_arrs=True)[channel].reshape(data[spatial_axis].shape)
    #compute combined background using region of overlap as a reference for magnitude
    overlap_magnitude = np.average(ph.roi(data, {0:spectral_reference_range, 1:spatial_reference_range}, return_arrs=True)[channel])
    background = spectral_background*spatial_background
    spectral_range = [__at(data[spectral_axis].points, spectral_reference_range[0]), __at(data[spectral_axis].points, spectral_reference_range[1])]
    spatial_range = [__at(data[spatial_axis].points, spatial_reference_range[0]), __at(data[spatial_axis].points, spatial_reference_range[1])]
    overlap_background = np.average(background[spectral_range[0]:spectral_range[1],spatial_range[0]:spatial_range[1]])
    overlap_ratio = overlap_magnitude/overlap_background
    background *= overlap_ratio
    #create background-subtracted channel
    if name is None:
        name = f"{channel}_bkgsub_combined"
    data.create_channel(name, values=data[channel][:]-background, units=data[channel].units)
    if data[channel].signed:
        data[name].signed = True
    if create_background_channel:
        data.create_channel(f"bkg_combined_{channel}", values=background, units=data[channel.units])
        if data[channel].signed:
            data[name].signed = True
    
    if talkback:
        print(f"subtracted combined background from channel {channel} of data {data.natural_name}")

def plot_profile(data, channel, **kwargs):
    
    #convert axis/channel indices to natural names
    channel, = ph.parse_args(data, channel, dtype='Channel')

    #set parameters for plotting from kwargs
    params = {
        'ROI' : None,
        'xticks' : None,
        'yticks' : None,
        'cbar_ticks' : None,
        'xlabel' : None,
        'ylabel' : None,
        'cbar_label' : None,
        "contrast" : None,
        "vrange" : None,
        "reference_lines" : None,
        "title" : None
        }
    params.update(styles.profile)
    
    if data[channel].signed:
        params["cmap"] = mpl.cm.RdBu_r   
    params.update(**kwargs)
    
    if params["ROI"] is not None:
        out = ph.roi(data, params["ROI"])
    else:
        out = data
    
    #determine range to be plotted
    if params["vrange"] is None:
        if params["contrast"] is None:
            vrange = ph.get_range(out, reference_key=channel)
        else:
            vrange = ph.contrast(out, channel, params["contrast"])
    else:
        vrange = params["vrange"]

    #setup x axis
    if params["xlabel"] is None:
        try:
            params["xlabel"] = out.variables[0].attrs['label']
        except KeyError:
            params["xlabel"] = 'spectrum'

    #setup y axis
    if params["ylabel"] is None:
        try:
            params["ylabel"] = out.variables[1].attrs['label']
        except KeyError:
            params["ylabel"] = 'y'

    #setup colorbar label
    if params["cbar_label"] is None:
        try:
            params["cbar_label"] = out[channel].attrs['label']
        except KeyError:
            params["cbar_label"] = 'signal'

    #setup plot frame
    fig, ax = plt.subplots(figsize=(params['fig_width'], params['fig_height']))
    
    #plot data
    xgrid, ygrid = np.meshgrid(out.axes[0][:], out.axes[1][:])
    #array needs to be transposed before passing to pcolormesh because apparently no matplotlib devs thought about what arrays look like
    try:
        mesh = ax.pcolormesh(xgrid, ygrid, np.transpose(out[channel][:]), cmap=params["cmap"], vmin=vrange[0], vmax=vrange[1])
    except TypeError:
        mesh = ax.pcolormesh(xgrid, ygrid, out[channel][:], cmap=params["cmap"], vmin=vrange[0], vmax=vrange[1])
    
    ax.set_xlabel(params["xlabel"])
    ax.set_ylabel(params["ylabel"])
    
    if params["title"] is not None:
        ax.set_title(params["title"])
    
    if params["reference_lines"] is not None:
        if type(params["reference_lines"]) is not list:
            params["reference_lines"] = [params["reference_lines"]]
        for reference_line in params["reference_lines"]:
            ax.axvline(x=reference_line, linewidth=1, color='grey', linestyle='--', alpha=0.25)
    
    #set ticks
    if not params["xticks"] and params['xticks'] is not None:
        ax.set_xticks([])
        
    if not params["yticks"] and params['yticks'] is not None:
        ax.set_yticks([])

    if not params['cbar_ticks'] and params['cbar_ticks'] is not None:
        ticks = []
    elif params['cbar_ticks'] is None:
        ticks = np.linspace(vrange[0], vrange[1], num=11)
    else:
        ticks = params['cbar_ticks']
    
    # plot colorbar
    cbar = plt.colorbar(mesh)
    cbar.set_ticks(ticks)
    cbar.set_label(params["cbar_label"])
    
def plot_decomposition(data, non_spatial_axis, spatial_axis, channel, **kwargs):
    #convert axis/channel indices to natural names
    non_spatial_axis, spatial_axis = ph.parse_args(data, non_spatial_axis, spatial_axis)
    channel, = ph.parse_args(data, channel, dtype='Channel')

    #set parameters for plotting from kwargs
    params = {
        "ROI" : None,
        "binning" : None,
        "xrange" : None,
        "vrange" : None,
        "yscale" : 'linear',
        "reference_lines" : None,
        "xticks" : True,
        "yticks" : True,
        "title" : None
        }
    params.update(styles.decomposition)
    params.update(**kwargs)

    #extract ROI
    if params["ROI"] is not None:
        out = ph.roi(data, params["ROI"])
    else:
        out = data
      
    #identify spatial ranges for indexing
    xrange = ph.get_range(out, reference_key=non_spatial_axis, dtype='Axis')

    #setup plot frame
    fig, ax = plt.subplots(figsize=(params['fig_width'], params['fig_height']))
    
    #plot data
    if params["binning"] is not None:
        if params["binning"] == 'average':
            arr_out = np.sum(out[channel][:], axis=1)/np.count_nonzero(np.sum(out[channel][:], axis=0))
        if params["binning"] == 'sum':
            arr_out = np.sum(out[channel][:], axis=1)
        #determine range to be plotted
        vrange = ph.vrange(arr_out, out[channel].signed, window=1)

        ax.plot(out[non_spatial_axis].points, arr_out,
            params["marker"], linewidth=params["linewidth"], alpha=1, color=params["color"])
    else:
        #determine range to be plotted
        vrange = ph.get_range(out, reference_key=channel)
        for i in range(out[spatial_axis].size):
            if np.sum(out[channel][:,i]) != 0:
                ax.plot(out[non_spatial_axis][:].flatten(), out[channel][:,i], 
                        params["marker"], linewidth=params["linewidth"], alpha=params["alpha"], color=params["color"])
                
    if out[channel].signed:
        ax.axhline(y=0, color='black', linewidth=1)
                
    if params["reference_lines"] is not None:
        if type(params["reference_lines"]) is not list:
            params["reference_lines"] = [params["reference_lines"]]
        for reference_line in params["reference_lines"]:
            ax.axvline(x=reference_line, zorder=0, linewidth=1, color='grey', linestyle='--', alpha=0.5)
    
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
    
    # label plot
    ax.set_xlabel(params["xlabel"])
    ax.set_ylabel(params["ylabel"])
    if params["title"] is not None:
        ax.set_title(params["title"])
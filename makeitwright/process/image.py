__name__ = "image"
__author__ = "Chris Roy, Song Jin Research Group, Dept. of Chemistry, University of Wisconsin - Madison"

#import
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
import processhelpers as proc
import styles 

def get_pixel_location(data, pixel):
    """
    Get the axis coordinates of an exact pixel in an image.

    Arguments
    ---------
    data : WrightTools Data - The image.
    pixel : tuple (x,y) - The pixel coordinates.

    Returns
    -------
    tuple (x,y) - The location of the pixel in axis coordinates.
    """
    return (data.axes[0].points[pixel[0]], data.axes[1].points[pixel[1]])

def remove_background(data, channel, threshold=0.5, negative=False, return_mask=False, max_ref_count=500):
    channel, = proc.parse_args(data, channel, dtype='Channel')
        
    ch_arr = data[channel][:]
    if max_ref_count > ch_arr.size:
        max_ref_count = int(ch_arr.size/10)
    ordered = np.sort(ch_arr.flatten())[-max_ref_count:]
    ch_max = np.average(ordered)
    bkg = np.where(ch_arr < threshold*ch_max, 0, ch_arr)
    mask = np.where(bkg > 0, 1, bkg)
    if negative:
        mask = 1-mask
    nobkg = ch_arr * mask

    data.create_channel(name=channel+"_nobkg", values=nobkg, units=data[channel].units)
    data[channel+'_nobkg'].signed = data[channel].signed
    if return_mask:
        data.create_channel(name=data[channel].natural_name+"_mask", values=mask)

def plot_image(data, channel, **kwargs):
    #convert axis/channel indices to natural names
    channel, = proc.parse_args(data, channel, dtype='Channel') 

    #set parameters for plotting from kwargs
    params = {
        "ROI" : None,
        "vrange" : None,
        "contrast" : None,
        "crosshairs" : None,
        "xticks" : None,
        "yticks" : None,
        "title" : None
        }
    params.update(styles.image)
    if data[channel].signed:
        params["cmap"] = mpl.cm.RdBu_r
    params.update(**kwargs)
    
    #extract ROI
    if params["ROI"] is not None:
        out = proc.roi(data, params["ROI"])
    else:
        out = data

    if params["xlabel"] is None:
        try:
            params["xlabel"] = out.variables[0].attrs['label']
        except KeyError:
            params["xlabel"] = 'x'

    if params["ylabel"] is None:
        try:
            params["ylabel"] = out.variables[1].attrs['label']
        except KeyError:
            params["ylabel"] = 'y'

    #determine range to be plotted
    if params["vrange"] is None:
        if params["contrast"] is None:
            vrange = proc.get_range(out, reference_key=channel)
        else:
            vrange = proc.contrast(out, channel, params["contrast"])
    else:
        vrange = params["vrange"]
    
    #setup plot frame
    fig, ax = plt.subplots(figsize=(params['fig_width'], params['fig_height']))
    
    #plot data
    xgrid, ygrid = np.meshgrid(out.axes[0][:], out.axes[1][:])
    ax.pcolormesh(xgrid, ygrid, np.transpose(out[channel][:]), cmap=params["cmap"], vmin=vrange[0], vmax=vrange[1])
    
    if params["crosshairs"] is not None:
        h, v = params["crosshairs"]
        if h is not None:
            ax.axvline(x=h, linewidth=1, color='white', linestyle='--', alpha=0.5)
        if v is not None:
            ax.axhline(y=v, linewidth=1, color='white', linestyle='--', alpha=0.5)

    if not params["xticks"] and params['xticks'] is not None:
        ax.set_xticks([])
    if not params["yticks"] and params['yticks'] is not None:
        ax.set_yticks([])

    ax.set_xlabel(params["xlabel"])
    ax.set_ylabel(params["ylabel"])

    ax.set_aspect("equal")

    if params["title"] is not None:
        ax.set_title(params["title"])
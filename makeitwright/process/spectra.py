__name__ = 'spectra'

import numpy as np
import WrightTools as wt
from matplotlib import pyplot as plt
import processhelpers as proc
import styles

"""
Data objects with one axis alone
"""

def plot_spectra(data, **kwargs):
    if type(data) is wt.Collection:
        data = [data[key] for key in data]
    if type(data) is not list:
        data = [data]
        
    #set parameters for plotting from kwargs
    params = {
        "plot_type" : "line",
        "xscale" : "linear",
        "xticks" : True,
        "yscale" : "linear",
        "yticks" : True,
        "axis" : 0,
        "channel" : -1,
        "ROI" : None,
        "xrange" : None,
        "vrange" : None,
        "offset" : 0,
        "reference_lines" : None,
        "title" : None,
        "background_color" : 'default'
        }
    params.update(styles.spectra)
    params.update(**kwargs)
        
    signed=False
    
    #parse color parameters to plot
    if type(params["colors"]) is list:
        colors = params["colors"]
        if len(params["colors"]) < len(data):
            q, r = divmod(len(data), len(colors))
            colors = q*colors+colors[:r]
    else: 
        try:
            colors = params["colors"](np.linspace(0,1,len(data)))
        except:
            colors = [params["colors"] for i in range(len(data))]
    
    #setup plot frame
    fig, ax = plt.subplots(figsize=(params['fig_width'], params['fig_height']))
    
    for i in range(len(data)):
        #convert axis/channel indices to natural names
        axis, = proc.parse_args(data[i], params["axis"])
        
        if params['channel']=='prompt': #kill this if else if all your code suddenly stops working
            channel, = proc.parse_args(data[i], input(f'select channel from {data[i].natural_name}: {[ch.natural_name  for ch in data[i].channels]} '),  dtype='Channel')
        else:    
            channel, = proc.parse_args(data[i], params["channel"],  dtype='Channel')
        if data[i][channel].signed:
            signed=True
    
        #extract ROI
        if params["ROI"] is not None:
            out = proc.roi(data[i], params["ROI"])
        else:
            out = data[i]
    
        #plot data
        if params["plot_type"] == "line":
            ax.plot(out[axis][:],out[channel][:]+i*params["offset"], 
                    linewidth=params["linewidth"], alpha=params["alpha"], color=colors[i])
        if params["plot_type"] == "scatter":
            ax.scatter(out[axis][:],out[channel][:]+i*params["offset"],
                    marker=params["marker"], alpha=params["alpha"], color=colors[i], s=params["marker_size"])
    
    if signed:
        ax.axhline(y=0, color='black', linewidth=1)
    
    if params["reference_lines"] is not None:
        if type(params["reference_lines"]) is not list:
            params["reference_lines"] = [params["reference_lines"]]
        for reference_line in params["reference_lines"]:
            ax.axvline(x=reference_line, zorder=0, linewidth=1, color='grey', linestyle='--', alpha=0.5)
    
    #adjust plot frame
    if params["xrange"] is not None:
        xrange = params["xrange"]
    else:
        xrange = proc.get_range(*data, reference_key=params["axis"], dtype='Axis')
    if params["xscale"] == 'log' and xrange[0]<=0:
        xrange[0] = 0.001
    ax.set_xlim(*xrange)
    ax.set_xscale(params["xscale"])
    if params["xlabel"] is None:
        try:
            params["xlabel"] = out[axis].attrs['label']
        except KeyError:
            params["xlabel"] = 'x'
    ax.set_xlabel(params["xlabel"])
    if not params["xticks"]:
        ax.set_xticks([])
    
    if params["vrange"] is not None:
        vrange = params["vrange"]
    else:
        vrange = proc.get_range(*data, reference_key=params["channel"], offset=params["offset"])
    if params["yscale"] == 'log' and vrange[0]<=0:
        vrange[0] = 0.01
    ax.set_ylim(*vrange)
    ax.set_yscale(params["yscale"])
    if params["ylabel"] is None:
        try:
            params["ylabel"] = out[channel].attrs['label']
        except KeyError:
            params["ylabel"] = 'y'
    ax.set_ylabel(params["ylabel"])
    if not params["yticks"]:
        ax.set_yticks([])

    if params["background_color"] != 'default':
        if params["background_color"] == 'transparent' or params["background_color"] is None:
            ax.set_alpha(0)
        else:
            ax.set_facecolor(params["background_color"])
        fig.set_alpha(0)

    if params["title"] is not None:
        ax.set_title(params["title"])

    plt.show()
        
def plot_tandem(d1,d2, figsize=(2.6,1), axis=0, channels=(-1,-1), 
                xticks=True, yticks=[True,True], xlabel="wavelength (nm)", ylabels=["reflectance","absorbance"], 
                xrange=[400,650], vranges=[(0,1),(0,1)], colors=['coral','royalblue'], 
                linewidth=1, reference_lines=None):
    #setup plot frame
    fig, ax1 = plt.subplots(figsize=figsize)
    ax2 = ax1.twinx()

    #convert axis/channel indices to natural names
    axis1, = proc.parse_args(d1, axis)
    axis2, = proc.parse_args(d2, axis)
    channel1, = proc.parse_args(d1, channels[0],  dtype='Channel')
    channel2, = proc.parse_args(d2, channels[1],  dtype='Channel')
    
    #plot data
    ax1.plot(d1[axis1][:],d1[channel1][:], linewidth=linewidth, color=colors[0])
    ax2.plot(d2[axis1][:],d2[channel2][:], linewidth=linewidth, color=colors[1])
    
    if reference_lines is not None:
        if type(reference_lines) is not list:
            reference_lines = [reference_lines]
        for line in reference_lines:
            ax1.axvline(x=line, zorder=0, linewidth=1, color='grey', linestyle='--', alpha=0.5)
    
    #adjust plot frame
    if xrange is None:
        xrange = proc.get_range(*[d1,d2], reference_key=axis, dtype='Axis')
    ax1.set_xlim(*xrange)
    ax1.set_xlabel(xlabel)
    if not xticks:
        ax1.set_xticks([])
    
    for i, v in enumerate(vranges):
        if v is None:
            if i==0:
                vranges[i] = proc.get_range(d1, reference_key=channel1, offset=0)
            if i==1:
                vranges[i] = proc.get_range(d2, reference_key=channel2, offset=0)
    ax1.set_ylim(*vranges[0])
    ax2.set_ylim(*vranges[1])
    ax1.set_ylabel(ylabels[0])
    ax2.set_ylabel(ylabels[1])
    if not yticks[0]:
        ax1.set_yticks([])
    if not yticks[1]:
        ax2.set_yticks([])
import numpy as np
import WrightTools as wt
import makeitwright.core.styles as styles

from .core import spectralprofile, hyperspectral

def central_wavelength(data):
    raise NotImplementedError


def plot_image(data, channel, **kwargs):
    params = {}
    try:
        unit = data['wl'].units
    except KeyError:
        unit = data.constants[0].units
    if unit == 'wn':
        params.update(styles.image_horiba_Raman)
    else:
        params.update(styles.image_horiba_PL)
    params.update(**kwargs)
    
    if len(data.axes) == 3:
        hyperspectral.plot_image(data, channel, **params)
    else:
        spectralprofile.plot_image(data, channel, **params)


def plot_profile(data, channel, profile_axis='y', **kwargs):
    params = {}
    try:
        unit = data['wl'].units
    except KeyError:
        unit = data.constants[0].units

    if data.axes[1].natural_name == 't':
        params.update(styles.profile_horiba_timed_series)
    elif unit == 'wn':
        params.update(styles.profile_horiba_Raman)
    else:
        params.update(styles.profile_horiba_PL)
    params.update(**kwargs)
    
    if len(data.axes) == 3:
        hyperspectral.plot_profile(data, profile_axis, channel, **params)
    else:
        spectralprofile.plot_profile(data, channel, **params)
    

def plot_decomposition(data, channel, **kwargs):
    params = {}
    try:
        unit = data['wl'].units
    except KeyError:
        unit = data.constants[0].units
    if unit == 'wn':
        params.update(styles.decomposition_horiba_Raman)
    else:
        params.update(styles.decomposition_horiba_PL)
    params.update(**kwargs)
    
    if len(data.axes) == 3:
        hyperspectral.plot_decomposition(data, 0, 1, 2, channel, **params)
    else:
        spectralprofile.plot_decomposition(data, 0, 1, channel, **params)

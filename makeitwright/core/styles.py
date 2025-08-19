import matplotlib.cm as cm
import cmocean


beckerhickl_transient = {
    "xlabel" : "time (ns)",
    "vreflines" : 0,
    "marker" : '.',
    "markersize" : 3
    }

profile = {
    "fig_width" : 4,
    "fig_height" : 3,
    "cmap" : cm.viridis,
    "xlabel" : None,
    "ylabel" : None,
    "cbar_label" : None,
    "cbar_ticks" : None
    }

profile_andor = {
    "fig_width" : 4,
    "fig_height" : 3
    }

profile_horiba = {
    "fig_width" : 6.5,
    "cbar_ticks" : None
    }

profile_horiba_PL = profile_horiba|{
    "fig_height" : 6.5,
    "cmap" : cm.viridis,
    "xlabel" : "wavelength (nm)",
    "ylabel" : "y (µm)",
    "cbar_label" : "PL intensity (cps)"
    }

profile_horiba_Raman = profile_horiba|{
    "fig_height" : 6.5,
    "cmap" : cm.inferno,
    "xlabel" : "Raman shift (cm\u207b\u2071)",
    "ylabel" : "y (µm)",
    "cbar_label" : "scattering intensity (cps)"
    }

profile_horiba_timed_series = profile_horiba|{
    "fig_height" : 10,
    "cmap" : cm.viridis,
    "xlabel" : "wavelength (nm)",
    "ylabel" : "excitation time (s)",
    "cbar_label" : "PL intensity (cps)"
    }

profile_iontof = {
    "fig_width" : 6.5,
    "fig_height" : 3.5,
    "cmap" : cmocean.cm.matter,
    "xlabel" : "distance (µm)",
    "ylabel" : "sputtering time (s)",
    "cbar_label" : "SI counts",
    "cbar_ticks" : None
    }

image = {
    "fig_width" : 4,
    "fig_height" : 4,
    "cmap" : cm.Greys_r,
    "xlabel" : None,
    "ylabel" : None,
    "cbar_label" : "signal (a.u.)",
    "ticks" : 'auto'
    }

image_andor = {
    "fig_width" : 4,
    "fig_height" : 4,
    "xlabel" : None,
    "ylabel" : None
    }

image_horiba = {
    "fig_width" : 6.5,
    "fig_height" : 6.5,
    "ticks" : 'auto'
    }

image_horiba_PL = image_horiba|{
    "cmap" : cm.viridis,
    "xlabel" : "x (µm)",
    "ylabel" : "y (µm)",
    "cbar_label" : "PL intensity (cps)"
    }

image_horiba_Raman = image_horiba|{
    "cmap" : cm.inferno,
    "xlabel" : "x (µm)",
    "ylabel" : "y (µm)",
    "cbar_label" : "scattering intensity (cps)"
    }

image_iontof = {
    "fig_width" : 6.5,
    "fig_height" : 6.5,
    "cmap" : cmocean.cm.matter,
    "xlabel" : "x (µm)",
    "ylabel" : "y (µm)",
    "cbar_label" : "SI counts",
    "ticks" : 'auto'
    }

decomposition = {
    "fig_width" : 6.5,
    "fig_height" : 3.0,
    "linewidth" : 1.0,
    "marker" : '-',
    "color" : 'black',
    "alpha" : 0.01,
    "xlabel" : "spectrum (a.u.)",
    "ylabel" : "signal (a.u.)"
    }

decomposition_andor = {
    "fig_width" : 6.5,
    "fig_height" : 3.5,
    "linewidth" : 1.0,
    "marker" : '-',
    "color" : 'red',
    "alpha" : 0.01,
    "xlabel" : "wavelength (nm)",
    }

decomposition_andor_A = decomposition_andor|{
    "ylabel" : "absorbance"
    }

decomposition_andor_PL = decomposition_andor|{
    "ylabel" : "PL intensity (cps)"
    }

decomposition_andor_R = decomposition_andor|{
    "ylabel" : "reflectance"
    }

decomposition_andor_RR0 = decomposition_andor|{
    "ylabel" : "reflection contrast"
    }

decomposition_andor_T = decomposition_andor|{
    "ylabel" : "transmittance"
    }

decomposition_horiba = {
    "fig_width" : 6.5,
    "fig_height" : 3.5,
    "linewidth" : 1.0,
    "marker" : '-',
    "color" : 'red',
    "alpha" : 0.01
    }

decomposition_horiba_PL = decomposition_horiba|{
    "xlabel" : "wavelength (nm)",
    "ylabel" : "PL intensity (cps)",
    }

decomposition_horiba_Raman = decomposition_horiba|{
    "xlabel" : "Raman shift (cm\u207b\u2071)",
    "ylabel" : "scattering intensity (cps)",
    }

decomposition_iontof = {
    "fig_width" : 6.5,
    "fig_height" : 2.0,
    "linewidth" : 1.0,
    "marker" : '.',
    "color" : 'red',
    "alpha" : 1,
    "xlabel" : "sputtering time (s)",
    "ylabel" : "SI counts"
    }

spectra = {
    "plot_type" : "line",
    "fig_width" : 4,
    "fig_height" : 3,
    "linewidth" : 2,
    "marker" : '.',
    "alpha" : 1,
    "marker_size" : 5,
    "colors" : cm.Set1,
    "xlabel" : None,
    "ylabel" : None
    }

spectra_TRPL = spectra|{
    "plot_type" : "scatter",
    "yscale" : "log",
    "xlabel" : "t (ns)",
    "ylabel" : "norm. counts",
    "colors" : cm.Set2,
    "marker" : '.',
    "marker_size" : 3,
    "reference_lines" : 0
    }

spectra_XRD_pattern = spectra|{
    "fig_height" : 3,
    "marker" : 'o',
    "marker_size" : 3,
    "colors" : cm.Set1,
    "xlabel" : "diffraction angle (deg. 2\u03B8)",
    "ylabel" : "intensity (a.u.)"
    }
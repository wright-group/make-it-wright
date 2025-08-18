import pathlib
import matplotlib as mpl
import makeitwright as mw


roi = mw.helpers.roi
parse = mw.parsers.parse
andor = mw.andor
becker = mw.beckerhickl
plot = mw.spectra.plot_spectra

fp = pathlib.Path().expanduser().resolve() / r"Desktop/Research Data/Wright Table/Original"

# set plotting parameters
mpl.rcParams['font.sans-serif'] = "Arial"
mpl.rcParams['font.family'] = "sans-serif"
mpl.rcParams['font.size'] = 14
mpl.rcParams['figure.dpi'] = 300
mpl.rcParams['lines.linewidth'] = 4
mpl.rcParams['pcolor.shading'] = 'auto'
mpl.rcParams['figure.dpi'] = 150

if True: # Plot PL
    data = parse(fp, keywords='4 hr.asc')
    PL_ROI = roi(data, {'y': ([1021, 1047], 'average')})
    plot(PL_ROI, channel=0, xrange=[500, 850])
    PL_output = open(fp / '4hr.txt', 'w')
    PL_dataTrace = zip(PL_ROI.axes[0], PL_ROI.channels[0])
    for x in PL_dataTrace:
        PL_output.write(str(x[0])+'\t')
        PL_output.write(str(x[1])+'\n')
    PL_output.close()

if True: # Plot T/R/A
    data = parse(fp / 'For Chris/23_11_21/4ClPEASnI n1', keywords='Object 3')
    R = data[2]
    R_back = data[1]
    T = data[4]
    T_back = data[3]
    
    andor.compute_reflectance(R, R_back, dark_wavelength_range=None)
    y_profile = roi(R, {'wl': ([580, 750], 'sum')})                           # If need to check object area
    plot(y_profile)
    plot(R, channel=1, ROI={'y': ([1020, 1070], 'average')}, xrange=[580, 750]) #Currently at 10 x mag
    R_ROI = roi(R, {'y': ([1020, 1070], 'average')})
    R_output = open(fp / 'For Chris/23_11_21/4ClPEASnI n1/Object 3 R processed.txt', 'w')
    R_dataTrace = zip(R_ROI.axes[0], R_ROI.channels[1])
    for x in R_dataTrace:
        R_output.write(str(x[0])+'\t')
        R_output.write(str(x[1])+'\n')
    R_output.close()
    
    andor.compute_transmittance(T, T_back, dark_wavelength_range=None)
    # y_profile = roi(T, {'wl': ([400, 500], 'sum')})                           # If need to check object area
    # plot(y_profile)
    plot(T, channel=1, ROI={'y': ([1020, 1070], 'average')}, xrange=[580, 750])     # Current 10x mag, 100x mag 54-70
    T_ROI = roi(T, {'y': ([1020, 1070], 'average')})
    T_output = open(fp / 'For Chris/23_11_21/4ClPEASnI n1/Object 3 T Processed.txt', 'w')
    T_dataTrace = zip(T_ROI.axes[0], T_ROI.channels[1])
    for x in T_dataTrace:
        T_output.write(str(x[0])+'\t')
        T_output.write(str(x[1])+'\n')
    T_output.close()
    #
    andor.compute_absorbance(R, T)
    A_output = open(fp / 'For Chris/23_11_21/4ClPEASnI n1/Object 3 A processed.txt', 'w')
    A_ROI = roi(T, {'y': ([1020, 1070], 'average')})                        # A is channel 2 in both R and T data objects
    plot(R, channel=2, ROI={'y': ([1020, 1070], 'average')}, xrange=[580, 750]) #Current 10x mag. can add vrange
    A_dataTrace = zip(A_ROI.axes[0], A_ROI.channels[2])
    for x in A_dataTrace:
        A_output.write(str(x[0])+'\t')
        A_output.write(str(x[1])+'\n')
    A_output.close()

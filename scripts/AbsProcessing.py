
# Process Reflectance/Transmittance/Absorbance Data from Wright group

import pathlib
import makeitwright.process.andor as andor
from makeitwright.process.helpers import roi
from makeitwright.parsers import parse
from makeitwright.artists import setparams, setdpi
from makeitwright.spectra import plot_spectra as plot


setparams()
setdpi(150)

filepath = pathlib.Path().expanduser().resolve() / "Desktop/Research Data/Wright Table/Original/test"
filename_R = "PEAPbBr4 R"
filename_RBack = "PEAPbBr4 R Back"
filename_T = "PEAPbBr4 T"
filename_TBack = "PEAPbBr4 T Back"
obj = 10            # Objective magnification (5x, 10x, 50x, 100x)
ROI_lower = 1000    # Lower and upper bounds of ROI
ROI_upper = 1047
plotx_lower = 350   # Lower and upper bounds of built-in plot x-axis
plotx_upper = 500

# Read data
data_R = parse(filepath, objective=obj, keywords=filename_R + ".asc")
data_RBack = parse(filepath, objective=obj, keywords=filename_RBack + ".asc")
data_T = parse(filepath, objective=obj, keywords=filename_T + ".asc")
data_TBack = parse(filepath, objective=obj, keywords=filename_TBack + ".asc")

processType = input(f'\nEnter "T" to compute transmittance, "R" to compute reflectance, or "A" to compute absorbance: ')

if processType == 'R':
    andor.compute_reflectance(data_R, data_RBack, dark_wavelength_range=None)

    # Check object area
    areaCheck = input(f'Enter "1" to check the Region of Interest. Otherwise, press enter: ')
    if areaCheck == '1':
        y_profile = roi(data_R, {'wl': ([400, 800], 'sum')})
        plot(y_profile)
        con = input(f'Enter "1" to stop program. Otherwise, press enter to continue to reflectance processing: ')
        if con == '1':
            quit()

    plot(data_R, channel=1, ROI={'y': ([ROI_lower, ROI_upper], 'average')}, xrange=[plotx_lower, plotx_upper])  # Can add vrange=[ , ] (y axis scale)

    # Print data to text file
    R_ROI = roi(data_R, {'y': ([ROI_lower, ROI_upper], 'average')})
    R_output = open(filepath + '/' + filename_R + ' Processed.txt', 'w')
    R_dataTrace = zip(R_ROI.axes[0], R_ROI.channels[1])
    for x in R_dataTrace:
        R_output.write(str(x[0])+'\t')
        R_output.write(str(x[1])+'\n')
    R_output.close()

if processType == 'T':
    andor.compute_transmittance(data_T, data_TBack, dark_wavelength_range=None)

    # Check object area
    areaCheck = input(f'Enter "1" to check the Region of Interest. Otherwise, press enter: ')
    if areaCheck == '1':
        y_profile = roi(data_T, {'wl': ([400, 800], 'sum')})
        plot(y_profile)
        con = input(f'Enter "1" to stop program. Otherwise, press enter to continue to reflectance processing: ')
        if con == '1':
            quit()

    plot(data_T, channel=1, ROI={'y': ([ROI_lower, ROI_upper], 'average')}, xrange=[plotx_lower, plotx_upper])  # Can add vrange=[ , ] (y axis scale)

    # Print data to text file
    T_ROI = roi(data_T, {'y': ([ROI_lower, ROI_upper], 'average')})
    T_output = open(filepath + '/' + filename_T + ' Processed.txt', 'w')
    T_dataTrace = zip(T_ROI.axes[0], T_ROI.channels[1])
    for x in T_dataTrace:
        T_output.write(str(x[0])+'\t')
        T_output.write(str(x[1])+'\n')
    T_output.close()

if processType == 'A':
    andor.compute_reflectance(data_R, data_RBack, dark_wavelength_range=None)
    andor.compute_transmittance(data_T, data_TBack, dark_wavelength_range=None)
    andor.compute_absorbance(data_R, data_T)

    # Check object area
    areaCheck = input(f'Enter "1" to check the Region of Interest. Otherwise, press enter: ')
    if areaCheck == '1':
        y_profile = roi(data_R, {'wl': ([400, 800], 'sum')})
        plot(y_profile)
        con = input(f'Enter "1" to stop program. Otherwise, press enter to continue to reflectance processing: ')
        if con == '1':
            quit()

    plot(data_R, channel=2, ROI={'y': ([ROI_lower, ROI_upper], 'average')}, xrange=[plotx_lower, plotx_upper])  # Can add vrange=[ , ] (y axis scale)

    # Print data to text file
    A_ROI = roi(data_T, {'y': ([ROI_lower, ROI_upper], 'average')})
    A_output = open(filepath + '/' + filename_T + ' Abs Processed.txt', 'w')    #TO DO finish naming
    A_dataTrace = zip(A_ROI.axes[0], A_ROI.channels[2])
    for x in A_dataTrace:
        A_output.write(str(x[0])+'\t')
        A_output.write(str(x[1])+'\n')
    A_output.close()


# Process PL Data from Wright group

import pathlib
import makeitwright as mw


andor = mw.andor
roi = mw.helpers.roi
parse = mw.parsers.parse
plot = mw.spectra.plot_spectra


filepath = pathlib.Path().expanduser() / "Desktop" / "Research Data" / "Wright Table" / "Original" / "test"
filename = "PEAPbI on FPEASnI PL 77K 4 2 hour wait for cool"
obj = 10            # Objective magnification (5x, 10x, 50x, 100x)
ROI_lower = 1000    # Lower and upper bounds of ROI
ROI_upper = 1047
plotx_lower = 500   # Lower and upper bounds of built-in plot x-axis
plotx_upper = 800


# Read data
data = parse(filepath, objective=obj, keywords=filename + ".asc")


# Check object area
areaCheck = input(f'Enter "1" to check the Region of Interest. Otherwise, press enter: ')
if areaCheck == '1':
    andor.correct_PL_background(data, ybkg=[0, 20])
    y_profile = roi(data, {'wl': ([400, 800], 'sum')})
    plot(y_profile)
    con = input(f'Enter "1" to stop program. Otherwise, press enter to continue to PL processing: ')
    if con == '1':
        quit()


# Process PL data
PL_ROI = roi(data, {'y': ([ROI_lower, ROI_upper], 'average')})
plot(PL_ROI, channel=0, xrange=[plotx_lower, plotx_upper])      # Can add vrange=[ , ] (y axis scale)


# Print data to text file
writeFile = input(f'Enter "1" to output results to text file. Otherwise, press enter: ')
if writeFile == '1':
    PL_output = open(filepath + '/' + filename + ' Processed.txt', 'w')
    PL_dataTrace = zip(PL_ROI.axes[0], PL_ROI.channels[0])
    for x in PL_dataTrace:
        PL_output.write(str(x[0])+'\t')
        PL_output.write(str(x[1])+'\n')
    PL_output.close()

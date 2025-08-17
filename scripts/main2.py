import numpy as np
import matplotlib.cm as cms
import cmocean.cm as cmo
from scipy.signal import savgol_filter as savgol
from scipy.signal import medfilt2d
from scipy.optimize import curve_fit
from scipy.stats import pearsonr, spearmanr, ttest_ind
import WrightTools as wt

import andor
import beckerhickl
import spectralprofile

from processhelpers import show, roi, norm, set_label
from parsers import parse
from artists import setparams, setdpi
from spectra import plot_spectra as plot

setparams()
setdpi(150)

fp = ""  # filepath name to folder



#
# Plot PL
data = parse('C:/Users/kmfor/Desktop/Research Data/Wright Table/Original', keywords='4 hr.asc')
#andor.correct_PL_background(data, ybkg=[0, 20])
#y_profile = roi(data, {'wl': ([400, 800], 'sum')})                           # If need to check object area
#plot(y_profile)
PL_ROI = roi(data, {'y': ([1021, 1047], 'average')})
plot(PL_ROI, channel=0, xrange=[500, 850])
PL_output = open('C:/Users/kmfor/Desktop/Research Data/Wright Table/Original/4hr.txt', 'w')
PL_dataTrace = zip(PL_ROI.axes[0], PL_ROI.channels[0])
for x in PL_dataTrace:
    PL_output.write(str(x[0])+'\t')
    PL_output.write(str(x[1])+'\n')
PL_output.close()

#
# Plot T/R/A
# data = parse('C:/Users/kmfor/Desktop/Research Data/Wright Table/Original/For Chris/23_11_21/4ClPEASnI n1', keywords='Object 3')
# R = data[2]
# R_back = data[1]
# T = data[4]
# T_back = data[3]

#
# andor.compute_reflectance(R, R_back, dark_wavelength_range=None)
# y_profile = roi(R, {'wl': ([580, 750], 'sum')})                           # If need to check object area
#plot(y_profile)
#plot(R, channel=1, ROI={'y': ([1020, 1070], 'average')}, xrange=[580, 750]) #Currently at 10 x mag
# R_ROI = roi(R, {'y': ([1020, 1070], 'average')})
# R_output = open('C:/Users/kmfor/Desktop/Research Data/Wright Table/Original/For Chris/23_11_21/4ClPEASnI n1/Object 3 R processed.txt', 'w')
# R_dataTrace = zip(R_ROI.axes[0], R_ROI.channels[1])
# for x in R_dataTrace:
#     R_output.write(str(x[0])+'\t')
#     R_output.write(str(x[1])+'\n')
# R_output.close()
#
# andor.compute_transmittance(T, T_back, dark_wavelength_range=None)
# # y_profile = roi(T, {'wl': ([400, 500], 'sum')})                           # If need to check object area
# # plot(y_profile)
#plot(T, channel=1, ROI={'y': ([1020, 1070], 'average')}, xrange=[580, 750])     # Current 10x mag, 100x mag 54-70
# T_ROI = roi(T, {'y': ([1020, 1070], 'average')})
# T_output = open('C:/Users/kmfor/Desktop/Research Data/Wright Table/Original/For Chris/23_11_21/4ClPEASnI n1/Object 3 T Processed.txt', 'w')
# T_dataTrace = zip(T_ROI.axes[0], T_ROI.channels[1])
# for x in T_dataTrace:
#     T_output.write(str(x[0])+'\t')
#     T_output.write(str(x[1])+'\n')
# T_output.close()
# #
# andor.compute_absorbance(R, T)
# A_output = open('C:/Users/kmfor/Desktop/Research Data/Wright Table/Original/For Chris/23_11_21/4ClPEASnI n1/Object 3 A processed.txt', 'w')
# A_ROI = roi(T, {'y': ([1020, 1070], 'average')})                        # A is channel 2 in both R and T data objects
# plot(R, channel=2, ROI={'y': ([1020, 1070], 'average')}, xrange=[580, 750]) #Current 10x mag. can add vrange
# A_dataTrace = zip(A_ROI.axes[0], A_ROI.channels[2])
# for x in A_dataTrace:
#     A_output.write(str(x[0])+'\t')
#     A_output.write(str(x[1])+'\n')
# A_output.close()

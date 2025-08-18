import matplotlib as mpl
from pathlib import Path

import makeitwright as mw

andor = mw.andor
roi = mw.helpers.roi
parse = mw.parsers.parse
plot = mw.artists.plot


user_path = Path().expanduser().resolve()

# set plotting parameters
mpl.rcParams['font.sans-serif'] = "Arial"
mpl.rcParams['font.family'] = "sans-serif"
mpl.rcParams['font.size'] = 14
mpl.rcParams['figure.dpi'] = 150
mpl.rcParams['lines.linewidth'] = 4
mpl.rcParams['pcolor.shading'] = 'auto'


# Plot PL
# data = parse('C:/Users/Kristel/Desktop/Research Data/Wright Table/Original/4_3_23', keywords=' PL')
# print(data)
# data[0].print_tree()
# # y_profile = roi(data, {'wm': ([400, 500], 'sum')})                           # If need to check object area
# # plot(y_profile)
# data_trace = andor.get_PL(data, ybkg=[0, 100], yspot = [1070, 1090])           # 1080 common center
# plot(data_trace)
# # #xrange=[390, 540]

# Plot T/R/A
data = parse('C:/Users/Kristel/Desktop/Research Data/Wright Table/Original/4_3_23')
R = data[2]
R_back = data[1]
T = data[4]
T_back = data[3]

A = andor.microabsorbance(T, R, T_back, R_back)
T_ROI = roi(T, {'yindex': ([1070, 1090], 'average')})
R_ROI = roi(R, {'yindex': ([1070, 1090], 'average')})
A_ROI = roi(A, {'yindex': ([1070, 1090], 'average')})
# R.print_tree()
#plot(R_ROI, xrange = [390, 530])
# R_dataTrace = zip(R.axes[0], R.channels[1])
# for x in R_dataTrace:
#     print(str(x[0])+'\t')
#     print(str(x[1])+'\n')

T_output = open(user_path / 'Desktop/Research Data/Wright Table/Original/4_3_23/PEA-Cs-Pb-Br n=2 4_3 T processed.txt', 'w')
T_dataTrace = zip(T_ROI.axes[0], T_ROI.channels[1])
for x in T_dataTrace:
    T_output.write(str(x[0])+'\t')
    T_output.write(str(x[1])+'\n')
T_output.close()

R_output = open(user_path / 'Desktop/Research Data/Wright Table/Original/4_3_23/PEA-Cs-Pb-Br n=2 4_3 R processed.txt', 'w')
R_dataTrace = zip(R_ROI.axes[0], R_ROI.channels[1])
for x in R_dataTrace:
    R_output.write(str(x[0])+'\t')
    R_output.write(str(x[1])+'\n')
R_output.close()


A_output = open(user_path / 'Desktop/Research Data/Wright Table/Original/4_3_23/PEA-Cs-Pb-Br n=2 4_3 A processed.txt', 'w')
A_dataTrace = zip(A_ROI.axes[0], A_ROI.channels[0])
for x in A_dataTrace:
    A_output.write(str(x[0])+'\t')
    A_output.write(str(x[1])+'\n')
A_output.close()

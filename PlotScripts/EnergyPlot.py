#!/usr/bin/python3
#!/home/fano.icmmg/berendeev_e_a/anaconda3/bin/python3
#!/mnt/storage/home/eaberendeev/anaconda3/bin/python3
#!/mnt/storage/home/aaefimova/anaconda3/bin/python3

import time
import multiprocessing
import numpy as np
import os
import matplotlib
#matplotlib.use('Qt4Agg')
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.colors as col
import matplotlib.ticker as ticker
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import gridspec
from matplotlib import rc
import pprint
import collections


#чтение параметров расчёта
from ReadDataLib import *
from LibPlot import *

def SubPlot(x,y,fig):
	return fig.add_subplot(gs[y,x])


#корневой каталог для расчёта (где файл параметров)
WorkDir='../'

#словарь с системными параметрами
SystemParameters=ReadParameters()
energy = ReadEnergyFile(WorkDir,SystemParameters)


fig = plt.figure(figsize=(12, 4))

#на основе того, сколько было сортов частиц определить сетку для графиков
gs =gridspec.GridSpec(1,2) #первый аргумент - вертикальное количество, второй - горизонтальное (по горизонтале число сортов частиц + 2 столбца для полей

ax = SubPlot(0,0,fig)		
ax.set_title('Energy')
ax.set_xlabel('$time, 1/w_p$')
ax.set_ylabel('Energy')
formatter = ticker.ScalarFormatter ()
formatter.set_powerlimits((-3, 3))
ax.yaxis.set_major_formatter (formatter)

for sort in SystemParameters['Particles'].keys():
	ax.plot(energy["time"],energy["Area_"+sort], label =sort)
ax.plot(energy["time"],energy["Area_E^2"],label ="E field")
ax.plot(energy["time"],energy["Area_B^2"]-energy["Area_B^2"][0],label ="B field")
ax.legend(loc='best')

ax = SubPlot(1,0,fig)		

ax.set_title('Radiation')
ax.set_xlabel('$time, 1/w_p$')
ax.set_ylabel('Radiation')
ax.yaxis.set_major_formatter (formatter)

ax.plot(energy["time"],energy["Rad_VacZ1"],label ="Vacuum Z_left")
ax.plot(energy["time"],energy["Rad_VacZ2"],label ="Vacuum Z_left")
ax.plot(energy["time"],energy["Rad_R"],label ="R")
ax.legend(loc='best')
gs.tight_layout(fig)
plt.savefig('../Anime/Energy.png', format='png', dpi=150)
plt.close(fig)
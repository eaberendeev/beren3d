#!/usr/bin/python3
#!/home/fano.icmmg/berendeev_e_a/anaconda3/bin/python3
#!/mnt/storage/home/eaberendeev/anaconda3/bin/python3
#!/mnt/storage/home/aaefimova/anaconda3/bin/python3

#from mpi4py import MPI

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

phasecdict = {'red':(
(0.0, 1, 1),#white
(0.15, 0, 0),#0096ff
(0.35, 0, 0),#blue
(0.55, 1, 1),#ff9600
(0.75, 1, 1),
(1.0, 0.5859375, 0.5859375)),

'green':   (
(0.0, 1, 1),#white
(0.15, 0.5859375, 0.5859375),#0096ff
(0.35, 0, 0),#blue
(0.55, 0.5859375, 0.5859375),#ff9600
(0.75, 0, 0),
(1.0, 0, 0)),

'blue':    (
(0.0, 1, 1),#white
(0.15, 1, 1),#0096ff
(0.35, 1, 1),#blue
(0.55, 0, 0),#ff9600
(0.75, 0.5859375, 0.5859375),
(1.0, 0, 0))}



#comm = MPI.COMM_WORLD
#rank = comm.Get_rank()
#ProcNum=comm.Get_size()
#print('ProcNum=',ProcNum)
##gs = gridspec.GridSpec(3, 4,height_ratios=[1,0.3,1],width_ratios=[1, 0.05,0.5,0.5])

#корневой каталог для расчёта (где файл параметров)
WorkDir='../'

FieldsAmp=0.025#амплитуда палитры для 3D полей
FieldsAmp2D=0.1#амплитуда палитры для 2D полей
#WindowSize=1000#окно для усреднения эффективности излучения

#словарь с системными параметрами
SystemParameters=ReadParameters()

def SubPlot(x,y,fig,gs):
	return fig.add_subplot(gs[y,x])

print(SystemParameters)

def split_rad_top(data_coord,data_rad):
	for t in range(1000):
		coord = []
		rad = []
		for p in range(100):
	#		filename = WorkDir + '/Fields//Diag1D//PRadTop_' + '{03}'.format(t) + '_' + '{03}'.format(p)
			t_str = "%03d" %t #'{03}'.format(t)
			p_str = "%03d" %p #'{03}'.format(p)
			filename = WorkDir + 'Fields//Diag1D//PRadAvgTop_' + t_str + '_' + p_str #'{03}'.format(t) #+ '_' + '{03}'.format(p)
			try:
				f = open(filename, 'r')
			except OSError:
				continue
			      #print("Все файлы обработаны!\n")
			titles = f.readline()
			for line in f:
				tmp = line.split()
				coord.append(float(tmp[0])) 
				rad.append(float(tmp[1]))
			f.close()
		if len(rad) > 0:
			data_coord.append(coord)
			data_rad.append(rad)

	for t in range(len(data_rad)):
		t_str = "%03d" %t
		filename = WorkDir + 'Fields//Diag1D//PRadAvgTop_' + t_str #'{03}'.format(t) #+ '_' + '{03}'.format(p)
		f = open(filename,"w")
		f.write(titles+'\n')
		for x in range(len(data_coord[t])):
			f.write(str((data_coord[t][x])) + " " + str(data_rad[t][x]))
			f.write('\n')
		f.close

def read_rad(data_coord,data_rad,title):
	for t in range(1000):
		t_str = "%03d" %t
		filename = WorkDir + 'Fields//Diag1D//PRadAvg'+title+'_' + t_str 
		try:
			f = open(filename, 'r')
		except OSError:
			continue
		coord = []
		rad = []
		titles = f.readline()
		for line in f:
			tmp = line.split()
			if len(tmp) == 2:
				coord.append(float(tmp[0])) 
				rad.append(float(tmp[1]))
		f.close()
		if len(rad) > 0:
			data_coord.append(coord)
			data_rad.append(rad)

def plot_rad1D(data_coord,data_rad,axis,title):
	maxT = []
	for t in range(len(data_rad)):
		maxT.append(max(data_rad[t]))
	
	maxRad = max(maxT)

	for t in range(len(data_rad)):
		t_str = "%03d" %t
		fig = plt.figure(figsize=(10, 6))
		gs = gridspec.GridSpec(1,1) #первый аргумент - вертикальное количество, второй - горизонтальное (по горизонтале число сортов частиц + 2 столбца для полей

		ax = SubPlot(0,0,fig,gs)
		formatter = ticker.ScalarFormatter ()
		formatter.set_powerlimits((-3, 3))
		ax.yaxis.set_major_formatter (formatter)
		ax.set_ylim(0, maxRad)
		ax.set_xlabel(axis+'$\ c/w_p$')
		ax.set_ylabel('$Prad$')
		ax.plot(data_coord[t],data_rad[t],label = "from "+title)

		#im = ax.imshow(data_rad, cmap=cm,origin='lower',extent=[data_coord[0][0],data_coord[0][len(data_coord[0])-1],0,len(data_rad)],aspect='auto')

		ax.set_title("Radiation " + title)

		plt.savefig('../Anime/Radiation/rad1D_' + title + t_str +'.png', format='png', dpi=150)
		plt.close(fig)



def plot_rad2D(data_coord,data_rad,axis,title):
	fig = plt.figure(figsize=(10, 5))
	gs = gridspec.GridSpec(1,1) #первый аргумент - вертикальное количество, второй - горизонтальное (по горизонтале число сортов частиц + 2 столбца для полей

	ax = SubPlot(0,0,fig,gs)
	cm = col.LinearSegmentedColormap('phase',phasecdict,N=1024,gamma=1)
	im = ax.imshow(data_rad, cmap=cm,origin='lower',extent=[data_coord[0][0],data_coord[0][len(data_coord[0])-1],0,len(data_rad)],aspect='auto')

	ax.set_xlabel(axis+'$\ c/w_p$')
	ax.set_ylabel('$time,\ 1/w_p$')
	ax.set_title("Radiation " + title)
	divider = make_axes_locatable(ax)
	cax = divider.append_axes("bottom", size="8%", pad=0.6) #свойства colorbar'а
	cbar =fig.colorbar(im, cax=cax, orientation='horizontal')
	tick_locator = ticker.MaxNLocator(nbins=5)
	cbar.locator = tick_locator
	cbar.formatter.set_powerlimits((0, 0))
	locator=ticker.MaxNLocator(prune=None, nbins=5)
	cbar.ax.xaxis.set_ticks_position('bottom')
	ax.yaxis.set_major_locator(locator)
	locator=ticker.MaxNLocator(prune=None, nbins=5)
	ax.xaxis.set_major_locator(locator)
	cbar.update_ticks()

	plt.savefig('../Anime/Radiation/rad2D_' + title +'.png', format='png', dpi=150)
	plt.close(fig)


try:
	os.mkdir('../Anime/Radiation')
except OSError:
	print("Каталог для излучения существуют")

data_coord = []
data_rad = []

read_rad(data_coord,data_rad,"Top")
plot_rad2D(data_coord,data_rad,"z","Top")
plot_rad1D(data_coord,data_rad,"z","Top")




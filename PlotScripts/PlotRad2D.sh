#!/home/eberendeev/anaconda3/bin/python3
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
import re


#чтение параметров расчёта
from ReadDataLib import *
from LibPlot import *


#comm = MPI.COMM_WORLD
#rank = comm.Get_rank()
#ProcNum=comm.Get_size()
#print('ProcNum=',ProcNum)
##gs = gridspec.GridSpec(3, 4,height_ratios=[1,0.3,1],width_ratios=[1, 0.05,0.5,0.5])

#корневой каталог для расчёта (где файл параметров)
WorkDir='../'

StartTime = 150 #len(Files)-1 #-rank #StartNum+  rank
node = "201"
nodeX = "591"
#словарь с системными параметрами
SystemParameters=ReadParameters()

SystemParameters['2D_Diag_Fields']=['1','1','1','1','1','1'] #ReadParamFile(WorkDir)

def SubPlot(x,y,fig):
	return fig.add_subplot(gs[y,x])


tdelay = int(SystemParameters['TimeDelay2D']) * float( SystemParameters['Dt'])

try:
      os.mkdir('../Anime/2D')
except OSError:
      print("подкаталог для 2D диагностики существуют")



#print(i,Files)

print(SystemParameters)


while StartTime >=0:
	TimeStep=str(StartTime).zfill(3) 
	#имя файла результирующего
	GraphName='Rad'+TimeStep+'.png'
	#проверка на наличие уже построенного графика
	print('TimeStep=',TimeStep)

	fig = plt.figure(figsize=(16, 12))

	#на основе того, сколько было сортов частиц определить сетку для графиков
	gs =gridspec.GridSpec(3,3,height_ratios=[0.1,1,1]) #первый аргумент - вертикальное количество, второй - горизонтальное (по горизонтале число сортов частиц + 2 столбца для полей

			
	RadPlaneY = []
	try: 
		RadPlaneY.append(ReadRadFile2D(WorkDir,"PlaneY000_",SystemParameters,TimeStep))
	except BaseException:
		pass
	try:
		RadPlaneY.append(ReadRadFile2D(WorkDir,"PlaneY001_",SystemParameters,TimeStep))
	except BaseException:
		pass

	RadPlaneZ = [] #collections.OrderedDict()
	try: 
		RadPlaneZ.append(ReadRadFile2D(WorkDir,"PlaneZ000_",SystemParameters,TimeStep))
	except BaseException:
		pass
	try:
		RadPlaneZ.append(ReadRadFile2D(WorkDir,"PlaneZ001_",SystemParameters,TimeStep))
	except BaseException:
		pass
	RadCircle =  [] #collections.OrderedDict()
	try: 
		RadCircle.append(ReadRadFile2D(WorkDir,"Circle000_",SystemParameters,TimeStep))
	except BaseException:
		pass
	try:
		RadCircle.append(ReadRadFile2D(WorkDir,"Circle001_",SystemParameters,TimeStep))
	except BaseException:
		pass

		
	gs_y=1
	for data in RadPlaneZ:
		Plot2DRad(data,"xy","Plane Z ",SystemParameters,SubPlot(0,gs_y,fig),fig)
		gs_y+=1

	gs_y=1
	for data in RadPlaneY:
		Plot2DRad(data,"xz","Plane Z",SystemParameters,SubPlot(1,gs_y,fig),fig)
		gs_y+=1
	gs_y=1
	for data in RadCircle:
		Plot2DRad(data,"circle","Circle",SystemParameters,SubPlot(2,gs_y,fig),fig)
		gs_y+=1
		
	CurrentTime=int(TimeStep)*tdelay
	CurrentTimeDt=round(CurrentTime,1)

	plt.figtext(0.5, 0.96,  r'$t\cdot\omega_p = '+str(CurrentTimeDt)+'$',bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5),
					color='black',fontsize=18,ha='center')

	plt.tight_layout()
	plt.savefig('../Anime/2D/'+GraphName, format='png', dpi=150)
	plt.close(fig)

	StartTime = StartTime - 1 


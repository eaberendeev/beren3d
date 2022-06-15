#!/usr/bin/python3
#!/home/eberendeev/anaconda3/bin/python3
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

FieldsAmp=0.025 #амплитуда палитры для 3D полей
StartTime = 10 #len(Files)-1 #-rank #StartNum+  rank
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
	GraphName='Res'+TimeStep+'.png'
	#проверка на наличие уже построенного графика
	print('TimeStep=',TimeStep)

	fig = plt.figure(figsize=(int(SystemParameters['NumOfPartSpecies'])*4+8, 12))

	#на основе того, сколько было сортов частиц определить сетку для графиков
	gs =gridspec.GridSpec(4,int(SystemParameters['NumOfPartSpecies'])+3,height_ratios=[0.1,1,1,1]) #первый аргумент - вертикальное количество, второй - горизонтальное (по горизонтале число сортов частиц + 2 столбца для полей
	#gs1 =gridspec.GridSpec(2,1,height_ratios=[0.1,1]) #первый аргумент - вертикальное количество, второй - горизонтальное (по горизонтале число сортов частиц + 2 столбца для полей

			
	#считываем плотности
	DensDataXY= collections.OrderedDict()
	for sort in SystemParameters['Particles'].keys():
		DensDataXY.update(ReadDensFile2D(WorkDir,"PlaneZ",SystemParameters,sort,TimeStep))
	#считываем плотности
	DensDataXZ= collections.OrderedDict()
	for sort in SystemParameters['Particles'].keys():
		DensDataXZ.update(ReadDensFile2D(WorkDir,"PlaneY",SystemParameters,sort,TimeStep))
#
	PhaseData= collections.OrderedDict()
	for sort in SystemParameters['Particles'].keys():
		PhaseData.update(ReadPhaseFile(WorkDir,SystemParameters,sort,TimeStep))

		#PlotIntegMPI(WorkDir)
	FieldDataXY=ReadFieldsFile2D(WorkDir,"PlaneZ_"+node+"_","xy",SystemParameters,TimeStep)
	FieldDataXZ=ReadFieldsFile2D(WorkDir,"PlaneY_"+node+"_","xz",SystemParameters,TimeStep)
	FieldDataYZ=ReadFieldsFile2D(WorkDir,"PlaneX_"+nodeX+"_","yz",SystemParameters,TimeStep)
		#pprint.pprint(FieldData)
		#pprint.pprint(FieldData.keys())

		#2D плотности

		
	gs_x=0
	for sort in SystemParameters['Particles'].keys():
		Plot2Ddens(DensDataXZ[sort],"xz",sort,0,2*float(SystemParameters['Particles'][sort]['Dens']),sort,SystemParameters,SubPlot(gs_x,1,fig),fig)
		Plot2Ddens(DensDataXY[sort],"xy",sort,0,2*float(SystemParameters['Particles'][sort]['Dens']),sort,SystemParameters,SubPlot(gs_x,2,fig),fig)
		Plot2DPhase(PhaseData,sort,0,2*float(SystemParameters['Particles'][sort]['Dens']),sort,SystemParameters,SubPlot(gs_x,3,fig),fig)
		#Plot1Ddens(DensData,sort,'long',0,4*float(SystemParameters['Particles'][sort]['Dens']),'',SystemParameters,fig.add_subplot(gs[gs_y,gs_x], label="1"),fig)
		
		gs_x+=1
		#gs_x=0
		#for sort in SystemParameters['Particles'].keys():
		#	gs_x+=1
		#gs_x=0
		#gs_y=2
		#for sort in SystemParameters['Particles'].keys():
		#	gs_x+=1
		#gs_x=0
		#gs_y=3
#		for sort in SystemParameters['Particles'].keys():
#			Plot1Ddens(DensData,sort,'long',320,0,4*float(SystemParameters['Particles'][sort]['Dens']),'',SystemParameters,fig.add_subplot(gs[gs_y,gs_x], label="1"),fig)
#			gs_x+=1

		#если есть правый пучок, то фазовый портрет его vx совмещаем с портретом от левого пучка
		#if 'RIGHTBEAM' in SystemParameters['Particles'].keys():
	#		PhaseData['LEFTBEAM']['vx']=PhaseData['LEFTBEAM']['vx']+PhaseData['RIGHTBEAM']['vx']
		#gs_x=0#поля начинаются со второго справа столбика
		#gs_y=1#первая строка
	#	for sort in SystemParameters['Particles'].keys():
	#		if sort!='RIGHTBEAM':#
#		Plot2Dphase(PhaseData,'BeamLeft','vx',0,0,'BeamLeft',SystemParameters,SubPlot(gs_x,2,fig),fig)
	#		Plot2Dphase(PhaseData,sort,'vy',0,0,sort,SystemParameters,SubPlot(gs_x,3,fig),fig)
	#		gs_x+=1


	#	gs_x=0
	#	for sort in SystemParameters['Particles'].keys():
	#		Plot1DField(FieldData,'Ex',int(Ny*0.5),-FieldsAmp2D,FieldsAmp2D,'',SystemParameters,fig.add_subplot(gs[-1,gs_x], label="2", frame_on=False),fig)
	#		gs_x+=1
	
	#	gs_x=-3 #поля начинаются со второго справа столбика
	#	gs_y=1 #первая строка
	#	for sort in FieldData.keys():
	#		if sort!= 'Bx':
     #                           Plot2Dfields(FieldData,sort,-FieldsAmp,FieldsAmp,sort,SystemParameters,SubPlot(gs_x,gs_y,fig),fig)
     #                           gs_y+=1
	#		if gs_y==4:
	#			gs_x+=1
	#			gs_y=1
	#		if sort== 'Bx':
	#			Plot2Dfields(FieldData,sort,float(SystemParameters['UniformB'][0])-FieldsAmp,float(SystemParameters['UniformB'][0])+FieldsAmp,sort,SystemParameters,SubPlot(gs_x,gs_y,fig),fig)
	#			gs_y+=1

	Plot2Dfields(FieldDataXZ,"xz","Ex",-FieldsAmp,FieldsAmp,"Ex",SystemParameters,SubPlot(-3,1,fig),fig)
	Plot2Dfields(FieldDataXZ,"xz","Ey",-FieldsAmp,FieldsAmp,"Ey",SystemParameters,SubPlot(-2,1,fig),fig)
	Plot2Dfields(FieldDataXZ,"xz","Ez",-FieldsAmp,FieldsAmp,"Ez",SystemParameters,SubPlot(-1,1,fig),fig)
	#Plot2Dfields(FieldDataXY,"xy","Ex",-FieldsAmp,FieldsAmp,"Ex",SystemParameters,SubPlot(-3,2,fig),fig)
	#Plot2Dfields(FieldDataXY,"xy","Ey",-FieldsAmp,FieldsAmp,"Ey",SystemParameters,SubPlot(-2,2,fig),fig)
	#Plot2Dfields(FieldDataXY,"xy","Ez",-FieldsAmp,FieldsAmp,"Ez",SystemParameters,SubPlot(-1,2,fig),fig)
	Plot2Dfields(FieldDataYZ,"yz","Ex",-FieldsAmp,FieldsAmp,"Ex",SystemParameters,SubPlot(-3,2,fig),fig)
	Plot2Dfields(FieldDataYZ,"yz","Ey",-FieldsAmp,FieldsAmp,"Ey",SystemParameters,SubPlot(-2,2,fig),fig)
	Plot2Dfields(FieldDataYZ,"yz","Ez",-FieldsAmp,FieldsAmp,"Ez",SystemParameters,SubPlot(-1,2,fig),fig)

	Plot2Dfields(FieldDataXZ,"xz","Bx",float(SystemParameters['UniformB'][0])-FieldsAmp,float(SystemParameters['UniformB'][0])+FieldsAmp,"Bx",SystemParameters,SubPlot(-3,3,fig),fig)
	Plot2Dfields(FieldDataXZ,"xz","By",-FieldsAmp,FieldsAmp,"By",SystemParameters,SubPlot(-2,3,fig),fig)
	Plot2Dfields(FieldDataXZ,"xz","Bz",-FieldsAmp,FieldsAmp,"Bz",SystemParameters,SubPlot(-1,3,fig),fig)

	CurrentTime=int(TimeStep)*tdelay
	CurrentTimeDt=round(CurrentTime,1)

	plt.figtext(0.5, 0.96,  r'$t\cdot\omega_p = '+str(CurrentTimeDt)+'$',bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5),
					color='black',fontsize=18,ha='center')

	plt.tight_layout()
	plt.savefig('../Anime/2D/'+GraphName, format='png', dpi=150)
	plt.close(fig)

#	GraphName='Rad'+TimeStep+'.png'

#	fig = plt.figure(figsize=(8, 6))
#	Rad = ReadRadFile2D(WorkDir,"",SystemParameters,TimeStep)

                #fig1 = plt.figure(figsize=16, 10))
#	MaxTime=int(TimeStep)*tdelay
#	MaxTimeDt=round(MaxTime,1)

                #на основе того, сколько было сортов частиц определить сетку для графиков
#	gs =gridspec.GridSpec(2,1,height_ratios=[0.1,1]) #первый аргумент - вертикальное количество, второй - горизонтальное (по горизонтале число сортов частиц + 2 столбца для полей
#	plt.figtext(0.5, 0.96,  r'$t\cdot\omega_p = '+str(MaxTimeDt)+'$',bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5),
 #                                               color='black',fontsize=18,ha='center')

                #axEff.legend()
#	Plot2DRad(Rad,0,2,"Radiation",SystemParameters,SubPlot(0,1,fig),fig)

#	plt.tight_layout()
#	plt.savefig('../Anime/2D/'+GraphName, format='png', dpi=150)
#	plt.close(fig)


#	GraphName='RadPlane'+TimeStep+'.png'
#	fig = plt.figure(figsize=(16, 12))
#	Rad = ReadRadFile2D(WorkDir,"_xy1_",SystemParameters,TimeStep)
                #fig1 = plt.figure(figsize=16, 10))
#	MaxTime=int(TimeStep)*tdelay
#	MaxTimeDt=round(MaxTime,1)

                #на основе того, сколько было сортов частиц определить сетку для графиков
#	gs =gridspec.GridSpec(3,2,height_ratios=[0.1,1,1]) #первый аргумент - вертикальное количество, второй - горизонтальное (по горизонтале число сортов частиц + 2 столбца для полей
#	plt.figtext(0.5, 0.96,  r'$t\cdot\omega_p = '+str(MaxTimeDt)+'$',bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5),
 #                                               color='black',fontsize=18,ha='center')

                #axEff.legend()
#	Plot2DRad(Rad,0,2,"Radiation xy1",SystemParameters,SubPlot(0,1,fig),fig)
#	Rad = ReadRadFile2D(WorkDir,"_xy2_",SystemParameters,TimeStep)
#	Plot2DRad(Rad,0,2,"Radiation xy2",SystemParameters,SubPlot(1,1,fig),fig)

#	Rad = ReadRadFile2D(WorkDir,"_xz1_",SystemParameters,TimeStep)
#	Plot2DRad(Rad,0,2,"Radiation xz1",SystemParameters,SubPlot(0,2,fig),fig)
	
#	Rad = ReadRadFile2D(WorkDir,"_xz2_",SystemParameters,TimeStep)
#	Plot2DRad(Rad,0,2,"Radiation xz2",SystemParameters,SubPlot(1,2,fig),fig)
	
#	plt.tight_layout()
#	plt.savefig('../Anime/2D/'+GraphName, format='png', dpi=150)
#	plt.close(fig)

#	GraphName='Fields'+TimeStep+'.png'

#	fig = plt.figure(figsize=(16, 8))
#	FieldDataC=ReadFieldsFile2D(WorkDir,"circle_","xy",SystemParameters,TimeStep)

                #fig1 = plt.figure(figsize=16, 10))
#	MaxTime=int(TimeStep)*tdelay
#	MaxTimeDt=round(MaxTime,1)

                #на основе того, сколько было сортов частиц определить сетку для графиков
#	gs =gridspec.GridSpec(2,4,height_ratios=[0.1,1]) #первый аргумент - вертикальное количество, второй - горизонтальное (по горизонтале число сортов частиц + 2 столбца для полей
#	plt.figtext(0.5, 0.96,  r'$t\cdot\omega_p = '+str(MaxTimeDt)+'$',bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5),
#                                                color='black',fontsize=18,ha='center')

                #axEff.legend()
#	Plot2Dfields(FieldDataC,"xz","Ex",-FieldsAmp,FieldsAmp,"Ex",SystemParameters,SubPlot(0,1,fig),fig)
#	Plot2Dfields(FieldDataC,"xz","By",-FieldsAmp,FieldsAmp,"By",SystemParameters,SubPlot(1,1,fig),fig)
#	Plot2Dfields(FieldDataC,"xz","Bz",-FieldsAmp,FieldsAmp,"Bz",SystemParameters,SubPlot(2,1,fig),fig)
#	Plot2Dfields2(FieldDataC,"xz","Bz",-FieldsAmp,FieldsAmp,"Pointing",SystemParameters,SubPlot(3,1,fig),fig)

#	plt.tight_layout()
#	plt.savefig('../Anime/2D/'+GraphName, format='png', dpi=150)
#	plt.close(fig)



	StartTime = StartTime - 1 


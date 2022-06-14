#!/usr/bin/python3


import time
import multiprocessing
import numpy as np
import os
import matplotlib
matplotlib.use('Qt4Agg')
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


#корневой каталог для расчёта (где файл параметров)
WorkDir='../'

FieldsAmp=0.025#амплитуда палитры для 3D полей
FieldsAmp2D=0.25#амплитуда палитры для 2D полей
WindowSize=1000#окно для усреднения эффективности излучения


def SubPlot(x,y,fig):
	return fig.add_subplot(gs[y,x])

TimeStep='0522'	
TimeStep='Last'#если это выбрано, то будет построен последний насчитанный файл


#чтение системных параметров
SystemParameters=ReadParamFile(WorkDir)
#вычислить число знаков в именах файлов
FileNameSignNum=len(str(int(int(SystemParameters['Max_Time'][0])/int(SystemParameters['Diagn_Time'][0]))))

if TimeStep=='Last':
	files=os.listdir(WorkDir+'Fields/DiagXY/')#получили список файлов в каталоге с полями (они есть всегда)
	#print(files)
	TimeStep=files[-1][-FileNameSignNum:]
	

print('TimeStep=',TimeStep)


fig = plt.figure(figsize=(16*2, 9*2))

#на основе того, сколько было сортов частиц определить сетку для графиков
gs =gridspec.GridSpec(5,len(SystemParameters['Particles'].keys())+2,height_ratios=[0.1,1,1,1,1])#первый аргумент - вертикальное количество, второй - горизонтальное (по горизонтале число сортов частиц + 2 столбца для полей

pprint.pprint(SystemParameters['Particles'].keys())
#pprint.pprint(SystemParameters)
Nx=int(SystemParameters['Nx'])
Ny=int(SystemParameters['Ny'])
dt=float(SystemParameters['dt'][0])


#считываем плотности
DensData= collections.OrderedDict()
for sort in SystemParameters['Particles'].keys():
	DensData.update(ReadDensFile(WorkDir,SystemParameters,sort,TimeStep))

DensVmax={}
for sort in SystemParameters['Particles'].keys():
	if(SystemParameters['Particles'][sort]['DensType'][0]=='0'):
		DensVmax[sort]=2*float(SystemParameters['Particles'][sort]['DensType'][1]);
	if(SystemParameters['Particles'][sort]['DensType'][0]=='1'):
		DensVmax[sort]=1.5*float(SystemParameters['Particles'][sort]['DensType'][2]);
	pprint.pprint(SystemParameters['Particles'][sort]['DensType'])
	pprint.pprint(DensVmax)

PhaseData= collections.OrderedDict()
for sort in SystemParameters['Particles'].keys():
	PhaseData.update(ReadPhaseFile(WorkDir,SystemParameters,sort,TimeStep))

#чтение данных о полях
FieldData=ReadFieldsFile(WorkDir,SystemParameters,TimeStep)

pprint.pprint(FieldData.keys())

#3D плотности
gs_x=0
for sort in SystemParameters['Particles'].keys():
	Plot2Ddens(DensData,sort,0,DensVmax[sort],sort,SystemParameters,SubPlot(gs_x,1,fig),fig)
	gs_x+=1

#если есть правый пучок, то фазовый портрет его vx совмещаем с портретом от левого пучка
if 'RIGHTBEAM' in SystemParameters['Particles'].keys():
	PhaseData['LEFTBEAM']['vx']=PhaseData['LEFTBEAM']['vx']+PhaseData['RIGHTBEAM']['vx']
gs_x=0#поля начинаются со второго справа столбика
gs_y=1#первая строка
for sort in SystemParameters['Particles'].keys():
	if sort!='RIGHTBEAM':
		Plot2Dphase(PhaseData,sort,'vx',0,0,sort,SystemParameters,SubPlot(gs_x,2,fig),fig)
	Plot2Dphase(PhaseData,sort,'vy',0,0,sort,SystemParameters,SubPlot(gs_x,3,fig),fig)
	gs_x+=1


gs_x=0
for sort in SystemParameters['Particles'].keys():
	Plot1Ddens(DensData,sort,'long',int(Ny*0.5),0,DensVmax[sort],'',SystemParameters,fig.add_subplot(gs[-1,gs_x], label="1"),fig)
	Plot1DField(FieldData,'Ex',int(Ny*0.5),-FieldsAmp2D,FieldsAmp2D,'',SystemParameters,fig.add_subplot(gs[-1,gs_x], label="2", frame_on=False),fig)
	gs_x+=1


gs_x=-2#поля начинаются со второго справа столбика
gs_y=1#первая строка
for sort in FieldData.keys():
	if sort!= 'Bx':
		Plot2Dfields(FieldData,sort,-FieldsAmp,FieldsAmp,sort,SystemParameters,SubPlot(gs_x,gs_y,fig),fig)
		gs_y+=1
	if gs_y==4:
		gs_x+=1
		gs_y=1
	if sort== 'Bx':
		Plot2Dfields(FieldData,sort,float(SystemParameters['UniformB'][0])-FieldsAmp,float(SystemParameters['UniformB'][0])+FieldsAmp,sort,SystemParameters,SubPlot(gs_x,gs_y,fig),fig)
		gs_y+=1


MaxTime=int(TimeStep)*int(SystemParameters['Diagn_Time'][0])

#каждый GPU пишет свой файл с интегральной энергией. Энергия полей в них одинаковая, а частиц - разная
IntegEnergy=ReadAllIntegFiles(WorkDir+'/Integs',SystemParameters,MaxTime)
MaxTimeDt=round(MaxTime*dt,1)

#pprint.pprint(IntegEnergy)
axIntError=SubPlot(-2,4,fig)

#print('MaxTime=',MaxTime)
time=np.arange(0,len(IntegEnergy['Error']))*dt

axIntError.plot(time,IntegEnergy['Error'],label='Error')
axIntError.set_xlim(left=0,right=MaxTimeDt)
axIntError.set_xlabel('$t/\omega_p$')

IntegEnergy['Efficiency']=IntegEnergy['EfficiencyTE']+IntegEnergy['EfficiencyTM']
axEff=SubPlot(-1,4,fig)

def moving_average(a, n):
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n
if(MaxTime>WindowSize):
	AvTE=moving_average(IntegEnergy['EfficiencyTE'],WindowSize)
	AvTM=moving_average(IntegEnergy['EfficiencyTM'],WindowSize)
	a=np.zeros(WindowSize-1)
	AvTE=np.insert(AvTE,0,a)
	AvTM=np.insert(AvTM,0,a)
else:
	AvTE=IntegEnergy['EfficiencyTE']
	AvTM=IntegEnergy['EfficiencyTM']
AvE=moving_average(IntegEnergy['Error'],WindowSize)
AvE=np.insert(AvE,0,a)
axEff.plot(time,AvTE*100,label='TE, %')
axIntError.plot(time,AvE,label='AvE')
axEff.plot(time,AvTM*100,label='TM, %')
axEff.plot(time,(AvTM+AvTE)*100,label='TE+TM, %')
axEff.set_xlim(left=0,right=MaxTimeDt)
axEff.set_ylim(bottom=0)
axEff.set_xlabel('$t/\omega_p$')

plt.figtext(0.5, 0.96,  r'$t\cdot\omega_p = '+str(round(MaxTimeDt,1))+'$',bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5),
				color='black',fontsize=18,ha='center')

plt.legend()

plt.tight_layout()
plt.show()

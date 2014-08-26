#! /usr/bin/env ipython
################################################################
##                                                            ##
##     Create time sereis of Av T, KE and MOC for spinup      ##
##                      Helen Burns                           ##
################################################################


''' Please give a run name in order to file away the output in
     ~/Figures/Runname
    This script will then calculate some spinup diagnostics 
    and draw a picture for all files in
    that folder and save the figures in ~/Figures/Folder.'''    
###################################################################        
################################################################### 
#--Import modules--#
from scipy.io import netcdf
import numpy as np
from pylab import clf, plot, show, floor, ceil, imshow
import matplotlib.pyplot as plt
import os
import csv
import sys
import glob
#--Take terminal inputs--#
if len(sys.argv) < 2:
   print '''\n    OOPSIE
            usage: Go to the folder with netcdf V output
            (named *all.nc with a grid folder located ../ or
            within the folder (for working directories) \n'''
   sys.exit(1)
OP = sys.argv[1]
#--Set folder structure--#
x=os.getcwd()
lists=glob.glob(x+'/*all.nc')
#--Main For loop--#
#For every .nc file in the folder
#Read in netcdf variables
#Decide resolution
#Find grid info
#Calculate the streamfunction
#Draw a picture and save it
if not lists:
   print '''\n    OOPSIE:
            no *all.nc folders here!!
            usage: Go to the folder with netcdf V output
            (named *all.nc with a grid folder located 
            within the folder (for working directories) \n'''
   sys.exit(1)
file2=netcdf.netcdf_file('grid.nc','r')
Zp1=file2.variables['Zp1']
Zp=Zp1[:]*1
file2=netcdf.netcdf_file('grid.nc','r')
Z=file2.variables['Z']
Z=Z[:]
Y=file2.variables['Yp1']
Y=Y[:]
## RMOC ##
def RMOC(file):
    file2read = netcdf.NetCDFFile(file,'r')
    lvrho=file2read.variables["LaVH1TH"]
    lvrho=lvrho[:]
    Y=file2read.variables['Yp1']
    Y=Y[:]
    time=file2read.variables['T']
    ti=time[:]
    dx=Y[1]-Y[0] # Find Resolution
    VT=np.nansum(lvrho*dx,axis=3) #integrate Vdx along x
    VTfdz=np.cumsum(VT[:,::-1,:],axis=1) #sum up the water column
    psi=np.nanmean(VTfdz[:,::-1,:],axis=0)/10**6 #Time average and put into Sv 
    #Take the minimum as the strength of ABBW Cells
    Psimin=np.nanmin((psi))
    Psimax=np.nanmax((psi))
    start=int(np.divide(ti[0],(86400*360)))#Find run start and stop times
    end=int(np.divide(ti[-1],(86400*360)))
    T=ti[0]
    if not Psi:
       return
    return Psimin,Psimax,T
## MOC ##
def MOC(file):
       file2read = netcdf.NetCDFFile(file,'r')
       V=file2read.variables['VVEL']
       V=V[:]
       #Find timestep and resolution
       time=file2read.variables['T']
       ti=time[:]
       dx=Y[1]-Y[0]
       Vtave=np.mean(V,axis = 0)
       Vtave[Vtave==0]=np.nan
       Vzone=np.nansum(Vtave,axis = 2)*dx
       dz=Zp[0:len(Zp)-1]-Zp[1:len(Zp)]
       psi = np.zeros((len(Zp),len(Y)))
       for j in range(len(Y)):
           for k in range(0,len(Z)-1):
               psi[k,j] = psi[k-1,j] + dz[k]*Vzone[k,j]
       y =Y/1000
       Psi=np.nanmax(abs(psi/10**6)) #Convert to Sv
       T=ti[0]
       start=int(np.divide(ti[0],(86400*360)))
       end=int(np.divide(ti[-1],(86400*360)))
       return Psi,T

## Now for AV T ##
def AvT(file):
       file2read = netcdf.NetCDFFile(file,'r')
       Temp=file2read.variables['THETA']
       temp=Temp[:]
       I=file2read.variables['T']
       time=I[:]
       file2=netcdf.netcdf_file('grid.nc','r')
       lm=file2.variables['HFacC']
       lm=lm[:]
       Temptav=np.nanmean(temp,axis=0)
       tempav=np.mean(Temptav*lm) 
       T=time[0]
       return tempav,T 

Psits=[]
TS=[]
for file in lists:
       Psi,T=MOC(file)
       Psits.append(Psi)
       TS.append(T)
Psits=np.array(Psits)
TS=np.array(TS)
year=TS/(86400*360)
a=np.array((year,Psits))
b=a[:,np.argsort(a[0,:])]
plt.plot(b[0,:],b[1,:])
plt.title("Max Mean MOC timeseries for years "+str(int(b[0,0]))+"-"+str(int(b[0,-1])))
plt.xlabel('Years')
plt.ylabel('$\overline{\psi}$ (sv)')
x=( os.path.expanduser('~')+"/Figures/"+OP)
plt.savefig(x+'/MOCtimeseries.png')
plt.clf()

tempts=[]
TS=[]
for file in lists:
       tempav,T=AvT(file)
       tempts.append(tempav)
       TS.append(T)
tempts=np.array(tempts)
TS=np.array(TS)
year=TS/(86400*360)
a=np.array((year,tempts))
b=a[:,np.argsort(a[0,:])]
plt.plot(b[0,:],b[1,:])
plt.title("Average domaine temperature timeseries for years "+str(int(b[0,0]))+"-"+str(int(b[0,-1])))
plt.xlabel('Years')
plt.ylabel('Temp $^o$C')
x=( os.path.expanduser('~')+"/Figures/"+OP)
plt.savefig(x+"/AvTtimeseries.png")
plt.clf()
x=os.getcwd()
lists=glob.glob(x+'/*Psi.nc')
psits=[]
psitsmin=[]
ts=[]
b=[]
a=[]
for file in lists:
       Psimin,Psimax,T=RMOC(file)
       psits.append(Psimax)
       psitsmin.append(Psimin)
       ts.append(T)
psits=np.array(psits)
psitsmin=np.array(psitsmin)
ts=np.array(ts)
if len(psits) != len(ts):
   print "WTF!"
   sys.exit(1)
year=ts/(86400*360)
a=np.array((year,psits))
b=a[:,np.argsort(a[0,:])]
plt.plot(b[0,:],b[1,:])
plt.title("Max Mean RMOC timeseries for years "+str(int(b[0,0]))+"-"+str(int(b[0,-1])))
plt.xlabel('Years')
plt.ylabel('$\psi ^{*}$ (sv)')
x=( os.path.expanduser('~')+"/Figures/"+OP)
plt.savefig(x+"/RMOCTSmax.png")
plt.clf()
a=np.array((year,psitsmin))
b=a[:,np.argsort(a[0,:])]
plt.plot(b[0,:],b[1,:])
plt.title("Min Mean RMOC timeseries for years "+str(int(b[0,0]))+"-"+str(int(b[0,-1])))
plt.xlabel('Years')
plt.ylabel('$\psi ^{*}$ (sv)')
x=( os.path.expanduser('~')+"/Figures/"+OP)
plt.savefig(x+"/RMOCTSmin.png")


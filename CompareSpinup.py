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
from numba import autojit
sys.path.append(os.path.expanduser('~')+"/bin/")
#--Take terminal inputs--#
''' We have 8 runs:                                                                                  
    3,10,30,100,300,1000,3000,and Closed '''
tau=[3,10,30,100,300,1000,3000]

#--Set folder structure--#
x=os.getcwd()

#--Main For loop--#
#For every .nc file in the folder
#Read in netcdf variables
#Decide resolution
#Find grid info
#Calculate the streamfunction
#Draw a picture and save it

#Check what files exist                                                                             
''' Let me know what files are missing and only bother                                               
    if at least 3 exists'''
check=0
runs=[]
for i in range(len(tau)):
    flist=str(tau[i])+'daynokpp/grid.nc'
    if not os.path.exists(flist):
       print 'WARNING: '+flist+' does not exist! (skipping this tau...)'
       check+=0
    else:
       check+=1
       runs.append(i)
if check<3:
   print "\n Less than 3 points to be plotted - no point. I'm stopping this..."
   sys.exit(1)
Runs=np.array(runs)
#print str(tau[Runs[0]])+'daynokpp/grid.nc'
file2=netcdf.netcdf_file(str(tau[Runs[0]])+'daynokpp/grid.nc','r')
Zp1=file2.variables['Zp1']
Zp=Zp1[:]*1
Z=file2.variables['Z']
Z=Z[:]
Y=file2.variables['Yp1']
Y=Y[:]
dx=Y[1]-Y[0] # Find Resolution
def regridy(Variable):
    Vc=(Variable[:,:,0:-1,:]+Variable[:,:,1::,:])/2
    return Vc
def regridx(Variable):
    Vc=(Variable[:,:,:,0:-1]+Variable[:,:,:,1::])/2
    return Vc
numba_regridy = autojit()(regridy)
numba_regridy.func_name = "numba_regridy"
numba_regridx = autojit()(regridx)
numba_regridx.func_name = "numba_regridx"
## EKE ##                                                                                            
def EKEf(U,V,dx):
       Vc=numba_regridy(V)
       Uc=numba_regridx(U)
       k=0.5*(Uc**2+Vc**2)
       EKE=np.nanmean(k)
       return EKE
## RMOC ##
def RMOC(file):
    file2read = netcdf.NetCDFFile(file,'r')
    lvrho=file2read.variables["LaVH1TH"]
    lvrho=lvrho[:]
    VT=np.nansum(lvrho*dx,axis=3) #integrate Vdx along x
    VTfdz=np.cumsum(VT[:,::-1,:],axis=1) #sum up the water column
    psi=np.nanmean(VTfdz[:,::-1,:],axis=0)/10**6 #Time average and put into Sv 
    #Take the minimum as the strength of ABBW Cells
    Psimin=np.nanmin((psi))
    Psimax=np.nanmax((psi))
    if not Psi:
       return
    return Psimin,Psimax
## MOC ##
def MOC(V,Zp):
       Vtave=np.mean(V,axis = 0)
       Vtave[Vtave==0]=np.nan
       Vzone=np.nansum(Vtave,axis = 2)*dx
       dz=Zp[0:len(Zp)-1]-Zp[1:len(Zp)]
       # Got rid of for loop here (much quicker!!)
       psi2=np.apply_along_axis(np.multiply,0,Vzone,dz)
       psi=np.cumsum(-psi2[::-1,:],axis=0)
       npad = ((0,1), (0,0))
       # Pad with zeros at bottom
       psi = np.pad(psi, pad_width=npad, mode='constant', constant_values=0)
       Psi=np.nanmax(abs(psi/10**6)) #Convert to Sv
       return Psi

## Now for AV T ##
def AvT(temp):
       Temptav=np.nanmean(temp,axis=0)
       tempav=np.mean(Temptav) 
       return tempav 
print "\n I've defined everything and I'm ready to start the Spin up diagnostics"

Tempplt=plt.subplot(321)
MOCplt=plt.subplot(322)
RMOCplt=plt.subplot(323)
RMOCminplt=plt.subplot(324)
EKEplt=plt.subplot(325)
for i in range(len(Runs)):
#for i in range(1):
    lists=glob.glob(x+'/'+str(tau[Runs[i]])+'daynokpp/*all.nc')
    Psits=[]
    EKEts=[]
    psits=[]
    tempts=[]
    psitsmin=[]
    TS=[]
    for file in lists:
           file2read = netcdf.NetCDFFile(file,'r')
           V=file2read.variables['VVEL']
           V=V[:]
           U=file2read.variables['UVEL']
           U=U[:]
           Temp=file2read.variables['THETA']
           temp=Temp[:]
           V=V[:]*1
           U=U[:]*1
           #Find timestep and resolution
           time=file2read.variables['T']
           ti=time[:]
           T=ti[0]
           a=list(file)
           a[-6]='P'
           a[-5]='s'
           a[-4]='i'
           A=''.join(a)
           Psi=MOC(V,Zp)
           tempav=AvT(temp)
           Psimin,Psimax=RMOC(A)
           EKE=EKEf(U,V,dx)
           EKEts.append(EKE)
           psits.append(Psimax)
           psitsmin.append(Psimin)
           tempts.append(tempav)
           Psits.append(Psi)
           TS.append(T)
           #print '.'
    EKEts=np.array(EKEts)
    psits=np.array(psits)
    psitsmin=np.array(psitsmin)
    Psits=np.array(Psits)
    tempts=np.array(tempts)
    TS=np.array(TS)
    year=TS/(86400*360)
    a1=np.array((year,Psits))
    b1=a1[:,np.argsort(a1[0,:])]
    a2=np.array((year,tempts))
    b2=a2[:,np.argsort(a2[0,:])]
    a3=np.array((year,psits))
    b3=a3[:,np.argsort(a1[0,:])]
    a4=np.array((year,psitsmin))
    b4=a4[:,np.argsort(a4[0,:])]
    a5=np.array((year,EKEts))
    b5=a5[:,np.argsort(a5[0,:])]
    MOCplt.plot(b1[0,:],b1[1,:])
    Tempplt.plot(b2[0,:],b2[1,:])
    RMOCplt.plot(b3[0,:],b3[1,:])
    RMOCminplt.plot(b4[0,:],b4[1,:])
    EKEplt.plot(b5[0,:],b5[1,:])
    print '\n'+str(tau[Runs[i]])+'day done! Only '+str(len(Runs)-i-1)+' to go...'
MOCplt.set_title(" Average MOC Strength)")
MOCplt.set_xlabel(r'Time')
MOCplt.set_ylabel('MOC (sv)')
EKEplt.set_title(" Average KE)")
EKEplt.set_xlabel(r'Time')
EKEplt.set_ylabel('KE')
Tempplt.set_xlabel('Time')
Tempplt.set_title(" Average Domain Temp)")
Tempplt.set_ylabel('Temp')
RMOCplt.set_title(" Average RMOC max)")
RMOCplt.set_xlabel(r'Time')
RMOCplt.set_ylabel('MOC (sv)')
RMOCminplt.set_title(" Average RMOC min)")
RMOCminplt.set_xlabel(r'Time')
RMOCminplt.set_ylabel('MOC (sv)')
runnames=[]
for i in range(len(Runs)):
    runnames.append(tau[Runs[i]])
lgd=plt.legend(runnames,bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
#MOCplt.legend(runnames,loc='upper center', shadow=True)
#EKEplt.legend(runnames,loc='upper center', shadow=True)
#Tempplt.plt.legend(runnames,loc='upper center', shadow=True)
#RMOCplt.legend(runnames,loc='upper center', shadow=True)
#RMOCminplt.legend(runnames,loc='upper center', shadow=True)
x=( os.path.expanduser('~')+"/Figures")
if not os.path.exists(x):
         os.makedirs(x)
#MOCplt.savefig(x+'/MOCspinup.png')
#EKEplt.savefig(x+'/KEspinup.png')
plt.tight_layout()
plt.savefig(x+'/spinup.png',bbox_extra_artists=(lgd,), bbox_inches='tight')
#RMOCplt.savefig(x+'/RMOCmaxspinup.png')
#RMOCminplt.savefig(x+'/RMOCminspinup.png')






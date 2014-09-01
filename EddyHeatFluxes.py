#! /usr/bin/env ipython                                                                              
################################################################                                     
##                                                            ##                                     
##        Script to create V'T' time average variable         ##                                    
##                        Helen Burns                         ##                                     
################################################################                                    
from scipy.io import netcdf
import numpy as np
import matplotlib.pyplot as plt
import os
#import csv
import sys
from pylab import *
from numba import autojit
#import glob
###
OP = sys.argv[1]
Comp='Mobilis'
x=os.getcwd()
#lists=glob.glob(x+'/'+str(OP)+'/*all.nc')
###
#Read in the timeaveraged files
file2 = netcdf.NetCDFFile(x+'/'+str(OP)+"/Tav.nc",'r') 
Temp=file2.variables['THETA']
Yc=file2.variables['Y']
Z=file2.variables['diag_levels']
X=file2.variables['X']
V=file2.variables['VVEL']
Yp1=file2.variables['Yp1']
#Regrid V to center points
def regrid(Variable):
    Vc=(Variable[:,:,0:-1,:]+Variable[:,:,1::,:])/2
    return Vc
numba_regrid = autojit()(regrid)
numba_regrid.func_name = "numba_regrid"
VC=numba_regrid(V[:]*1)
# Remove excess time dimension
VTAVC=np.squeeze(VC)
TTAV=np.squeeze(Temp[:]*1)
#Area matrix
#x+'/'+str(OP)+grd="/noc/msm/scratch/students/hb1g13/"+Comp+"/"+OP+"/grid/"
#os.chdir(grd)
file2=netcdf.netcdf_file(x+'/'+str(OP)+'/grid.nc','r')
Zp1=file2.variables['Zp1']
Zp=Zp1[:]*1
dx=5000
dz=Zp[0:len(Zp)-1]-Zp[1:len(Zp)]
dA=dz*5000
VT=VTAVC*TTAV
VTdA=np.zeros(shape(VT))
for k in range(len(dz)):
    VTdA[k,:,:]=VT[k,:,:]*dA[k]
MeanHF=1030*3985*(np.sum(np.sum(VTdA,axis=0),axis=1))/10**15
#Read in the VTtav files
fileVT = netcdf.NetCDFFile(x+'/'+str(OP)+"/VTbar.nc",'r') 
VTbar=fileVT.variables['VT']
VTdA=np.zeros(shape(VT))
for k in range(len(dz)):
    VTdA[k,:,:]=VTbar[k,:,:]*dA[k]
TotalHF=1030*3985*(np.sum(np.sum(VTdA,axis=0),axis=1))/10**15
#Read in the VTtav files
fileVTp = netcdf.NetCDFFile(x+'/'+str(OP)+"/VTprimebar.nc",'r') 
VTbar=fileVTp.variables['VT']
for k in range(len(dz)):
    VTdA[k,:,:]=VTbar[k,:,:]*dA[k]
EddyHF=1030*3985*(np.sum(np.sum(VTdA,axis=0),axis=1))/10**15
E,=plt.plot(EddyHF)
M,=plt.plot(MeanHF,'r')
T,=plt.plot(TotalHF,'k')
plt.title("Meridional Heat Fluxes ("+OP+")")
plt.ylabel("Heat Flux (PW)")
plt.xlabel("Meridional Distance (km)")
lgd=legend([E,M,T],["Eddy $V'T'$","Mean $\overline{V}\,\,\overline{T}$","Total $\overline{VT}$"],bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
x=( os.path.expanduser('~')+"/Figures/"+Comp+"/"+OP)
if not os.path.exists(x):
          os.makedirs(x)
plt.savefig(x+"/EMTplt.png",bbox_extra_artists=(lgd,), bbox_inches='tight')
vzone=np.mean(VTAVC,axis=2)
tzone=np.mean(TTAV,axis=2)
def gyrecomps(Tav1,Zone1,Tav2,Zone2):
    T1star=np.zeros(shape(Tav1))
    T2star=np.zeros(shape(Tav1))
    T1sz=np.zeros(shape(Tav1))
    T2sz=np.zeros(shape(Tav1))
    for i in range(len(Tav1[1,1,:])):
        T1star[:,:,i]=Tav1[:,:,i]-Zone1
        T2star[:,:,i]=Tav2[:,:,i]-Zone2
        T1sz[:,:,i]=T1star[:,:,i]*Zone2
        T2sz[:,:,i]=T2star[:,:,i]*Zone1
    return T1star, T2star, T1sz, T2sz
numba_gyrecomps= autojit()(gyrecomps)
numba_gyrecomps.func_name = "numba_gyrecomps"
[T1star, T2star, T1sz, T2sz]=numba_gyrecomps(TTAV,tzone,VTAVC,vzone)
vtzone=vzone*tzone
vtstar=T1star*T2star
def xdA(Term,dA):
    VTdA=np.zeros(shape(Term))
    for k in range(len(dz)):
        VTdA[k,:,:]=Term[k,:,:]*dA[k]
    return VTdA
numba_xdA= autojit()(xdA)
numba_xdA.func_name = "numba_xdA"
T2=numba_xdA(T1sz,dA)
T3=numba_xdA(T2sz,dA)
T4=numba_xdA(vtstar,dA)
T1=np.zeros(shape(vtzone))
for k in range(len(dz)):
    T1[k,:]=vtzone[k,:]*dA[k]
G1=200*1030*3985*(np.sum(T1,axis=0))/10**15
G2=1030*3985*(np.sum(np.sum(T2,axis=0),axis=1))/10**15
G3=1030*3985*(np.sum(np.sum(T3,axis=0),axis=1))/10**15
G4=1030*3985*(np.sum(np.sum(T4,axis=0),axis=1))/10**15
GyreHF=G2+G3+G4
MeanHF=G1
E,=plt.plot(EddyHF)
M,=plt.plot(MeanHF,'r')
G,=plt.plot(GyreHF,'g')
T,=plt.plot(TotalHF,'k')
plt.title("Meridional Heat Fluxes ("+OP+")")
plt.ylabel("Heat Flux (PW)")
plt.xlabel("Meridional Distance (km)")
lgd=legend([E,M,G,T],["Eddy $V'T'$","Mean <$\overline{V}><\overline{T}>$","Gyre $V^*T^*$","Total $\overline{VT}$"],bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.savefig(x+'/'+OP+"/EMTGplt.png",bbox_extra_artists=(lgd,), bbox_inches='tight')
EddyHF0=1030*3985*(np.sum(np.sum(VTdA[0:3,:,:],axis=0),axis=1))/10**15
EddyHF1=1030*3985*(np.sum(np.sum(VTdA[0:5,:,:],axis=0),axis=1))/10**15
EddyHF2=1030*3985*(np.sum(np.sum(VTdA[0:10,:,:],axis=0),axis=1))/10**15
EddyHF3=1030*3985*(np.sum(np.sum(VTdA[0:15,:,:],axis=0),axis=1))/10**15
EddyHF4=1030*3985*(np.sum(np.sum(VTdA[0:20,:,:],axis=0),axis=1))/10**15
EddyHF0=1030*3985*(np.sum(np.sum(VTdA[0:3,:,:],axis=0),axis=1))/10**15
E0,=plt.plot(EddyHF0,'k')
E1,=plt.plot(EddyHF1)
E2,=plt.plot(EddyHF2,'r')
E3,=plt.plot(EddyHF3,'g')
E4,=plt.plot(EddyHF4,'c')
E5,=plt.plot(EddyHF,'m')
plt.title("Meridional Eddy Fluxes ("+OP+"s)")
plt.ylabel("Heat Flux (PW)")
plt.xlabel("Meridional Distance (km)")
lgd=legend([E0,E1,E2,E3,E4,E5],[str(Zp[3]),str(Zp[5]),str(Zp[10]),str(Zp[15]),str(Zp[20]),str(Zp[30])],bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.savefig(x+"/"+OP+"/Eddysplt.png",bbox_extra_artists=(lgd,), bbox_inches='tight')

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 30 20:38:55 2024

@author: isaias
"""




# nohup /opt/anaconda3_titan/bin/python /home/isaias.ramirez/GNSSETAS_EM2.py &
from mpmath import mp
import numpy as np
import pandas as pd
import scipy as sp
import geopy.distance
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import time
from multiprocessing import Pool, cpu_count, get_context, Process
from scipy import interpolate
import plotly.io as pio
import plotly.graph_objs as go
pio.renderers.default='browser'
import plotly.express as px
from scipy.optimize import minimize
from scipy.stats import multivariate_normal
import os
try:
    os.chdir('/home/isaias.ramirez/')
except:
    os.chdir('/Users/isaias/Desktop/GNSS')

try :
    os. mkdir("./figures")
except:
    pass

plt.rcParams['figure.dpi'] = 300

############# Define useful functions

M0=4.3
def fintegrated(i,alpha,d,Lat,Lon):
    cov=(d*mp.exp(alpha*(Datos["Magnitude"][i]-M0)))*np.diag((1,1))
    cum1=multivariate_normal.cdf((latmax,lonmax), mean=(Lat,Lon), cov=cov)
    cum2=multivariate_normal.cdf((latmin,lonmax), mean=(Lat,Lon), cov=cov)
    cum3=multivariate_normal.cdf((latmax,lonmin), mean=(Lat,Lon), cov=cov)
    cum4=multivariate_normal.cdf((latmin,lonmin), mean=(Lat,Lon), cov=cov)
    return cum1-cum2-cum3+cum4


def fintegrated2(i,alpha,d,Lat,Lon,Lim):
    cov=(d*mp.exp(alpha*(Datos["Magnitude"][i]-M0)))*np.diag((1,1))
    cum1=multivariate_normal.cdf((Lat+Lim,Lon+Lim), mean=(Lat,Lon), cov=cov)
    cum2=multivariate_normal.cdf((Lat+Lim,Lon-Lim), mean=(Lat,Lon), cov=cov)
    cum3=multivariate_normal.cdf((Lat-Lim,Lon+Lim), mean=(Lat,Lon), cov=cov)
    cum4=multivariate_normal.cdf((Lat-Lim,Lon-Lim), mean=(Lat,Lon), cov=cov)
    return cum1-cum2-cum3+cum4


def klam(M,A,alpha,M0):
    return A*mp.exp(alpha*(M-M0))

def fint(x,y,M, d,alpha,M0):
    return 1/(2*np.pi*d*mp.exp(alpha*(M-M0)))*mp.exp(-(x**2+y**2)/(2*d*np.exp(alpha*(M-M0))))

def g(t,p,c):
    return (p-1)*mp.power(c,p-1)*mp.power(t+c,-p)*(t>0)

def RecuperaBloque(Longitude,Latitude):
    for i in np.arange(mu0.shape[0]):
        for j in np.arange(mu0.shape[1]):
            # print(i,j)
            # print(len(Datos[(DomX[j]<Datos["Longitude"])*(DomX[j+1]>=Datos["Longitude"])*(DomY[i]<Datos["Latitude"])*(DomY[i+1]>=Datos["Latitude"])]))
            Eleccion=(DomX[j]<=Longitude)*(DomX[j+1]>=Longitude)*(DomY[i]<=Latitude)*(DomY[i+1]>=Latitude)
            if Eleccion==1:
                return (i,j)

def lam(i):
    xObs,yObs,MObs=Datos.loc[i][["Latitude","Longitude","Magnitude"]]
    t=Fechas[i]
    Bloque=RecuperaBloque(yObs,xObs) 
    muxy=mu0[Bloque]+mu1[Bloque]*(PII[i,1]!=0)+mu2[Bloque]*(PII[i,2]!=0)+mu3[Bloque]*(PII[i,3]!=0)#+mu4[Bloque]*(PII[i,4])
    nu=0
    nu=np.zeros(0)
    for j in np.arange(i):
        Lon,Lat,M2=Datos.loc[j][["Longitude","Latitude","Magnitude"]]
        # nu+= klam(M2,A,alpha,M0)*fint(xObs-Lat,yObs-Lon,Mag, d,alpha,M0)*g((t-Fechas[j]).total_seconds()/sectoday,p,c)
        nu=np.hstack((nu,klam(M2,A,alpha,M0)*fint(xObs-Lat,yObs-Lon,M2, d,alpha,M0)*g((t-Fechas[j]).total_seconds()/sectoday,p,c))) #Puede que haya error, cotegar si es MObs o M2
    return muxy,nu

def likelihhod(thetalike):
    A,alpha,c,p,d=np.copy(thetalike)
    
    global lik        
    def lik(j):
        value=0
        xObs=Datos["Latitude"].loc[j]
        yObs=Datos["Longitude"].loc[j]
        t= mindate+ timedelta(seconds=Diferencias[j])
        Mag=Datos["Magnitude"].loc[j]
        intg=(1-mp.power(c,(p-1))*(((maxdate-t).total_seconds())/sectoday+c)**(1-p))
        intf=fintegrated(j,alpha,d,xObs,yObs)
        for i in np.arange(j,len(Datos)-1)+1:    
            Lat=Datos["Latitude"].loc[i]
            Lon=Datos["Longitude"].loc[i]
            # Mag2=Datos["Magnitude"].loc[i]
            value+= PIJ[i,j]*mp.log(klam(Mag,A,alpha,M0)*fint(xObs-Lat,yObs-Lon,Mag, d,alpha,M0)*g((Fechas[i]-t).total_seconds()/sectoday,p,c))
        value=value-klam(Mag,A,alpha,M0)*intg*intf#- 10000*(fintegrated2(j,alpha,d,xObs,yObs,0.5)-1)**2 #-(intf-1)**2-(intg-1)**2
        return value
    paralell= get_context("fork").Pool(cpu_count())
    results = paralell.map(lik, range(len(Datos)))
    paralell.close()
    value=np.sum(results)
    print(thetalike,"\n", value)
    global thetaIteracion
    thetaIteracion=np.copy(thetalike)
    indmax=Datos[Datos["Magnitude"]==np.max(Datos["Magnitude"])].index[0]
    LatMax=Datos["Latitude"].loc[indmax]
    LonMax=Datos["Longitude"].loc[indmax]
    # return value#-100000*((1-mp.power(c,(p-1))*((10+c)**(1-p)))-1)**2-100*(klam(7.6,A,alpha,M0)-15)**2-10000*(fintegrated2(indmax,alpha,d,LatMax,LonMax,1)-1)**2
    return value-10_000*((1-mp.power(c,(p-1))*((10+c)**(1-p)))-1)**2-10_000*(fintegrated2(indmax,alpha,d,LatMax,LonMax,1)-1)**2#-1*(klam(7.6,A,alpha,M0)-15)**2

################ Hyperparameter
mindate=datetime(2000, 1,1)
sectoday=86400
anio=2012  #2009 #2012 #Fijar en 2017
maxdate=datetime(anio, 1,1)
################ Read data GNSS
file_path = './CAYA_2000-2021_GAMIT.dat'
# file_path = '/Users/isaias/Desktop/PINO_GAMIT.dat'
GNSSData=np.ones(9)
with open(file_path, 'r') as file:
    for line in file:
        GNSSData=np.vstack((GNSSData,np.fromstring((line.strip()).replace("  ", " "), sep=" ")))
        # print(np.fromstring((line.strip()).replace("  ", " "), sep=" ") )

GNSSData=GNSSData[1:,:]
GNSSData=GNSSData[GNSSData[:,0]<anio]
for i in range(len(GNSSData)):
  if i==0:
    GNSSDate=datetime(int(GNSSData[i,0]),int(GNSSData[i,1]),int(GNSSData[i,2]))
  else :
    GNSSDate=np.hstack((GNSSDate,datetime(int(GNSSData[i,0]),int(GNSSData[i,1]),int(GNSSData[i,2]))))

GNSSDif=np.zeros(len(GNSSData))
for i in range(len(GNSSDif)):
    GNSSDif[i]=(GNSSDate[i]-mindate).total_seconds()

################ Read data Epicenter
Datos=pd.read_csv('./hypo_MASUB.csv',delimiter=",")
Datos=Datos.loc[ Datos["Year"]<anio]
Datos=Datos.loc[ Datos["Magnitude"]>=M0]
Datos=Datos.reset_index()
for i in range(len(Datos)):
    if i==0:
        Fechas=datetime(Datos["Year"][i], Datos["Month"][i], Datos["Day"][i], Datos["Hour"][i], Datos["Minute"][i], int(Datos["Second"][i]), int(str(Datos["Second"][i])[-1])*100000)
    else :
        Fechas=np.hstack((Fechas,datetime(Datos["Year"][i], Datos["Month"][i], Datos["Day"][i], Datos["Hour"][i], Datos["Minute"][i], int(Datos["Second"][i]), int(str(Datos["Second"][i])[-1])*100000)))

Diferencias=np.zeros(len(Fechas))
for i in range(len(Diferencias)):
    Diferencias[i]=(Fechas[i]-mindate).total_seconds()

############### Exploring data
# plt.plot(GNSSDif,GNSSData[:,3])    
# plt.scatter(Diferencias,np.log(Datos["Magnitude"]))
    
# plt.plot(GNSSDif,GNSSData[:,3])    
# plt.scatter(Diferencias,np.cumsum(np.log(Datos["Magnitude"]))) 
S=1 #Extra size of the square
latmin=np.min(Datos["Latitude"])-S
latmax=np.max(Datos["Latitude"])+S
lonmin=np.min(Datos["Longitude"])-S
lonmax=np.max(Datos["Longitude"])+S    
################ Read data Epicenter
x=GNSSDif/sectoday
y=GNSSData[:,3] 
# y=sp.signal.detrend(y) #Detrend

#############################
plt.plot(GNSSDate,y)
plt.xlabel("Date")
plt.ylabel("Displacement (cm)")

#############################

tck = interpolate.splrep(x, y, s=0.03) #############################Tener cuidado si se cambia el soporte, 2017 era 0.5
xnew = x
ynew = interpolate.splev(xnew, tck, der=0)
plt.plot(x,y)
plt.plot(xnew,ynew)
yder = interpolate.splev(xnew, tck, der=1)   # or BSpline(*tck)(xnew, 1)
plt.plot(x,yder)
plt.scatter(xnew,ynew)
plt.scatter(xnew[yder<0],(ynew[yder<0]))
plt.scatter(xnew[yder<0],100*np.abs(yder[yder<0]))
cumsum=np.cumsum(yder<0)
plt.plot(cumsum)
DominioInt=np.zeros(0,dtype=int)
flag=2
flagN=-1
for i in range(len(cumsum)-1):
    if cumsum[i]==cumsum[i+1] and flag!=1:
        flag=1
    if flag!=0 and  (cumsum[i]==cumsum[i+1]-1):
        flag=0
    if flagN!=flag:
        DominioInt=np.append(DominioInt,i)
        flagN=flag

DominioInt=np.unique(DominioInt[:100])
DominioInt=DominioInt[1:]
plt.plot(GNSSDate,y,color="orange")
plt.plot(GNSSDate,ynew)


plt.plot(GNSSDate,y)
plt.plot(GNSSDate,ynew)
# plt.scatter(GNSSDate[yder<0],100*np.abs(yder[yder<0]),s=0.5)
for i in range(len(DominioInt)-1):
    plt.axvline(GNSSDate[DominioInt][i],color="firebrick")
plt.xlabel("Date")
plt.ylabel("Displacement (cm)")

###############################################################################Algorithm 1 Fox

PII=np.zeros((len(Fechas),len(DominioInt)//2+1))
PII[:,0]=1/(np.arange(len(Fechas))+1)
for i in range(len(DominioInt)//2):
    PII[(GNSSDate[DominioInt][2*i]<Fechas) * (GNSSDate[DominioInt][2*i+1]>Fechas),0]=PII[(GNSSDate[DominioInt][2*i]<Fechas) * (GNSSDate[DominioInt][2*i+1]>Fechas),0]/2
    PII[(GNSSDate[DominioInt][2*i]<Fechas) * (GNSSDate[DominioInt][2*i+1]>Fechas),(i+1)]=PII[(GNSSDate[DominioInt][2*i]<Fechas) * (GNSSDate[DominioInt][2*i+1]>Fechas),0]

###############################################################################

PIJ=np.zeros((len(Fechas),len(Fechas)))
for i in range(len(Fechas)):
    for j in np.arange(i):
        PIJ[i,j]=1/(i+1)


Trench=pd.read_csv('./trench_cocos.csv', delimiter=',')



npartX=20
npartY=20
DomX=np.min(Datos["Longitude"])+np.arange(npartX)/(npartX-1)*(np.max(Datos["Longitude"])-np.min(Datos["Longitude"]))
DomY=np.min(Datos["Latitude"])+np.arange(npartY)/(npartY-1)*(np.max(Datos["Latitude"])-np.min(Datos["Latitude"]))
M=Datos["Magnitude"]

plt.xlim(DomX[0],DomX[-1])
plt.ylim(DomY[0],DomY[-1])
for i in range(len(DomX)):
    plt.axvline(DomX[i],color="firebrick",alpha=0.5)
    plt.axhline(DomY[i],color="firebrick",alpha=0.5)
plt.plot(Trench["x"],Trench["y"],color="black")
plt.scatter(Datos["Longitude"],Datos["Latitude"],s=np.exp(5*(M-np.min(M))/(np.max(M)-np.min(M))+1),alpha=0.5)
plt.ylabel("Latitude")
plt.xlabel("Longitude")

DeltaX=DomX[1]-DomX[0]
DeltaY=DomY[1]-DomY[0]


###############################################################################
mu0=np.zeros((npartY-1,npartX-1))
mu1=np.zeros((npartY-1,npartX-1))
mu2=np.zeros((npartY-1,npartX-1))
mu3=np.zeros((npartY-1,npartX-1))
mu4=np.zeros((npartY-1,npartX-1))


T0=(maxdate-mindate).total_seconds()/sectoday

for i in range(len(DominioInt)//2):
    T0=np.hstack((T0,(GNSSDate[DominioInt][2*i+1]-GNSSDate[DominioInt][2*i]).total_seconds()/sectoday))


###############################################################################
thetaopt=(0.05, #antes 0.04
          1.5,
          0.02,
          1.5,
          0.004)

A,alpha,c,p,d=thetaopt

Const=[[10**(-5), 10**(0)],#A 
[0.5, 2],#alpha
[10**(-3), 10**(1)],#c
[1.1, 3],#p
[0+10**(-10), 10**(-2)]]#d antes era -1 con este cambio no necesito penalizar la esperanza 


###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################


for s in range(40):
    for i in np.arange(mu0.shape[0]):
        for j in np.arange(mu0.shape[1]):
            # print(i,j)
            # print(len(Datos[(DomX[j]<Datos["Longitude"])*(DomX[j+1]>=Datos["Longitude"])*(DomY[i]<Datos["Latitude"])*(DomY[i+1]>=Datos["Latitude"])]))
            Eleccion=PII[(DomX[j]<Datos["Longitude"])*(DomX[j+1]>=Datos["Longitude"])*(DomY[i]<Datos["Latitude"])*(DomY[i+1]>=Datos["Latitude"])]
            #print(Eleccion)
            mu0[i,j]=np.sum(Eleccion[:,0])/(DeltaX*DeltaY*T0[0])
            mu1[i,j]=np.sum(Eleccion[:,1])/(DeltaX*DeltaY*T0[1])
            mu2[i,j]=np.sum(Eleccion[:,2])/(DeltaX*DeltaY*T0[2])
            mu3[i,j]=np.sum(Eleccion[:,3])/(DeltaX*DeltaY*T0[3])
            #mu4[i,j]=np.sum(Eleccion[:,4])/(DeltaX*DeltaY*T0[4])  ####Descomentar si es hasta 2017
    MAX=np.max(np.hstack((mu0,mu1,mu2,mu3))) #Falta mu4    
    MIN=np.min(np.hstack((mu0,mu1,mu2,mu3))) #Falta mu4    
    # plt.imshow(mu0[:,::-1], extent=[DomX[0],DomX[1],DomY[0],DomY[1]],aspect=DeltaX/DeltaY)
    fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(14, 8))
    axs[0, 0].set_title('mu0')
    im1=axs[0, 0].imshow(mu0[::-1,:], extent=[DomX[0],DomX[-1],DomY[0],DomY[-1]],aspect=DeltaX/DeltaY,vmin=np.min(MIN,0),vmax=MAX)
    plt.colorbar(im1)
    axs[0,0].scatter(Datos["Longitude"],Datos["Latitude"],s=np.exp(5*(M-np.min(M))/(np.max(M)-np.min(M))+1),alpha=0.1,color="red")
    axs[0, 1].set_title('mu1')
    axs[0, 1].imshow(mu1[::-1,:], extent=[DomX[0],DomX[-1],DomY[0],DomY[-1]],aspect=DeltaX/DeltaY,vmin=0,vmax=MAX)    
    axs[1, 0].set_title('mu3')
    axs[1, 0].imshow(mu3[::-1,:], extent=[DomX[0],DomX[-1],DomY[0],DomY[-1]],aspect=DeltaX/DeltaY,vmin=0,vmax=MAX)
    if np.sum(mu4!=0):
        axs[1, 1].set_title('mu4')
        axs[1, 1].imshow(mu4[::-1,:], extent=[DomX[0],DomX[-1],DomY[0],DomY[-1]],aspect=DeltaX/DeltaY,vmin=0,vmax=MAX)
    else:
        axs[1, 1].set_title('$P_{ii}$')
        axs[1, 1].hist(np.apply_along_axis(sum, 1,PII))
    
    axs[0, 2].set_title('mu2')
    axs[0, 2].imshow(mu2[::-1,:], extent=[DomX[0],DomX[-1],DomY[0],DomY[-1]],aspect=DeltaX/DeltaY,vmin=0,vmax=MAX)
    axs[1, 2].set_title('PIJ')
    axs[1, 2].imshow((PIJ[:100,:100]))
    fig.show()
    fig.tight_layout()
    fig.savefig("./figures/"+str(s)+"a.jpg", dpi=250)
    fig2, axs = plt.subplots(nrows=1, ncols=1, figsize=(14, 8))
    im1=axs.imshow((PIJ[:100,:100]))
    plt.colorbar(im1)
    fig2.show()
    fig2.tight_layout()
    fig2.savefig("./figures/"+str(s)+"b.jpg", dpi=250)
    fig3, axs = plt.subplots(nrows=1, ncols=1, figsize=(14, 8))
    axs.scatter(np.arange(len(PIJ)),np.apply_along_axis(sum, 0,PIJ))
    axs.set_title(thetaopt)
    fig3.show()
    fig3.tight_layout()
    fig3.savefig("./figures/"+str(s)+"c.jpg", dpi=250)
    
    ###############################################################################
    Optimizacion=minimize(lambda theta: -likelihhod(theta),thetaopt, bounds=Const, method="L-BFGS-B",tol=1e-6 )#"Nelder-Mead"  "Powell" "L-BFGS-B"
    print(Optimizacion)
    A,alpha,c,p,d=Optimizacion["x"]
    # Optimo=thetaIteracion
    # A,alpha,c,p,d=Optimo
    # A,alpha,c,p,d=[0.00531568, 0.6095514,  0.02656009 ,1.99870522 ,0.00608718] 
    thetaopt=(A,alpha,c,p,d)
    for i in range(len(PII)):
        xObs,yObs,MObs=Datos.loc[i][["Latitude","Longitude","Magnitude"]]
        t=Fechas[i]
        Bloque=RecuperaBloque(yObs,xObs)  #######Aqui podria haber error
        muxy,nu=lam(i)
        lamObs=muxy+np.sum(nu)
        
        PII[i,0]=mu0[Bloque]/(lamObs)#*(PII[i,0])
        PII[i,1]=mu1[Bloque]/(lamObs)*(PII[i,1]!=0)
        PII[i,2]=mu2[Bloque]/(lamObs)*(PII[i,2]!=0)
        PII[i,3]=mu3[Bloque]/(lamObs)*(PII[i,3]!=0)
        #PII[i,4]=mu4[Bloque]/(lamObs)*(PII[i,4]!=0)*(PII[i,4])
        for j in np.arange(i):
            PIJ[i,j]=nu[j]/lamObs
    np.savetxt("./figures/PII.csv", PII, delimiter=",")
    np.savetxt("./figures/PIJ.csv", PIJ, delimiter=",")
    np.savetxt("./figures/mu0.csv", mu0, delimiter=",")
    np.savetxt("./figures/mu1.csv", mu1, delimiter=",")
    np.savetxt("./figures/mu2.csv", mu2, delimiter=",")
    np.savetxt("./figures/mu3.csv", mu3, delimiter=",")
    np.savetxt("./figures/thetaopt.csv", thetaopt, delimiter=",")


np.apply_along_axis(sum, 1,PIJ)+np.apply_along_axis(sum, 1,PII)


########No usar la aproximacion, falta restar en likelihood los parametros



fig2 = px.imshow(PIJ,labels=dict(x="Index", y="Index"))
fig2.show()
fig2.write_html("./figures/Genealogy"+".html")


# fig2 = px.imshow(PII)
# fig2 = px.imshow(np.hstack((PII,np.zeros((len(PII),len(PII))))))
# fig2.show()
#fig2.write_html("./figures/Genealogy"+".html")




# Extension='/Users/isaias/Desktop/figures/'
# PII=np.loadtxt(Extension+'PII.csv', delimiter=',')
# PIJ=np.loadtxt(Extension+'PIJ.csv', delimiter=',')
# mu0=np.loadtxt(Extension+'mu0.csv', delimiter=',')
# mu1=np.loadtxt(Extension+'mu1.csv', delimiter=',')
# mu2=np.loadtxt(Extension+'mu2.csv', delimiter=',')
# mu3=np.loadtxt(Extension+'mu3.csv', delimiter=',')
# thetaopt=np.loadtxt(Extension+'thetaopt.csv', delimiter=',')

# A,alpha,c,p,d=thetaopt


# suma=np.sum(mu0*T0[0]+mu1*T0[1]+mu2*T0[2]+mu3*T0[3])*(DeltaX*DeltaY)
# for i in Datos.index[1:]:
#     xObs,yObs,MObs=Datos.loc[i][["Latitude","Longitude","Magnitude"]]
#     t=Fechas[i]
#     intg=(1-mp.power(c,(p-1))*(((maxdate-t).total_seconds())/sectoday+c)**(1-p))
#     intf=fintegrated(j,alpha,d,xObs,yObs)
#     suma+=klam(MObs,A,alpha,M0)*intg*intf

# print(suma,len(PIJ))

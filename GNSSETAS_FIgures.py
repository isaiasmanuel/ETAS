#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 11:47:05 2024

@author: isaias
"""

#Run GNSSETAS_EM2 before

from mpmath import mp
import numpy as np
import pandas as pd
import scipy as sp
import geopy.distance
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import time
from multiprocessing import Pool, cpu_count, get_context, Process

    
################ Hyperparameter

mindate=datetime(2000, 1,1)
sectoday=86400
anio=2017
maxdate=datetime(anio, 1,1)
################ Read data GNSS
file_path = './CAYA_2000-2021_GAMIT.dat'
# file_path = '/Users/isaias/Desktop/PINO_GAMIT.dat'

 

GNSSData=np.ones(9)
with open(file_path, 'r') as file:
    for line in file:
        GNSSData=np.vstack((GNSSData,np.fromstring((line.strip()).replace("  ", " "), sep=" ")))
        print(np.fromstring((line.strip()).replace("  ", " "), sep=" ") )
        
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

Datos=pd.read_csv('/Users/isaias/Desktop/hypo_MASUB.csv',delimiter=",")
Datos=Datos.loc[ Datos["Year"]<anio]


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
    
    
    
    
import plotly.io as pio
import plotly.graph_objs as go

pio.renderers.default='browser'
import plotly.express as px

    
################ Read data Epicenter

CoordStation=(17.0485,-100.2672) #CAYA
# CoordStation=(16.3928,-98.1273) #PINO
StationDist=np.ones(len(Datos))

for i in range(len(Datos)):
    StationDist[i]=geopy.distance.geodesic((Datos["Latitude"][i],Datos["Longitude"][i]),CoordStation).km


r=100 #np.inf
Selected=StationDist<r

from sklearn.neighbors import KernelDensity
kde =KernelDensity(kernel="gaussian",bandwidth=10_000_000).fit(np.reshape(Diferencias[Selected],(len(Diferencias[Selected]),1)))    
# plt.plot(GNSSDif,GNSSData[:,3])        
# plt.plot(Diferencias[Selected],10**(7)*np.exp(kde.score_samples(np.reshape(Diferencias[Selected],(len(Diferencias[Selected]),1))))) 
    

fig = px.line(x=GNSSDif/sectoday, y=GNSSData[:,3] )
#fig.add_trace(go.Scatter(x=Diferencias, y=Datos["Magnitude"], mode='markers'))
fig.add_trace(go.Scatter(x=Diferencias[Selected]/sectoday,y=10**(7)*np.exp(kde.score_samples(np.reshape(Diferencias[Selected],(len(Diferencias[Selected]),1)))), marker_size=5*np.exp(Datos[Selected]["Magnitude"]-4), mode='lines+markers'))
fig.show()
        

plt.hist(Datos[Selected]["Magnitude"])
 
    
    










import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

x=GNSSDif/sectoday
y=GNSSData[:,3] 
# y=sp.signal.detrend(y) #Detrend
plt.plot(GNSSDate,y)
plt.ylabel("Desplazamiento [cm]")

tck = interpolate.splrep(x, y, s=0.05)

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
    #i=338
    if cumsum[i]==cumsum[i+1] and flag!=1:
        flag=1
                
    if flag!=0 and  (cumsum[i]==cumsum[i+1]-1):
        flag=0
        
    if flagN!=flag:
        
        DominioInt=np.append(DominioInt,i)
        flagN=flag


DominioInt=np.unique(DominioInt[:100])
DominioInt=DominioInt[1:]


plt.scatter(xnew,ynew)
plt.scatter(xnew[yder<0],100*np.abs(yder[yder<0]))

for i in range(len(DominioInt)):
    plt.axvline(x[DominioInt][i])

x[DominioInt][0::2]
x[DominioInt][1::2]



plt.scatter(xnew[yder<0],np.abs(yder[yder<0]))




    
#############################    
    
    

import geopandas as gpd

M=Datos["Magnitude"]

Mapa1 = './2001.geojson'
Mapa2 = './2006.geojson'
Mapa3 = './2009.geojson'
Mapa4 = './2014.geojson'
df1 = gpd.read_file(Mapa1)
df2 = gpd.read_file(Mapa2)
df3 = gpd.read_file(Mapa3)
df4 = gpd.read_file(Mapa4)




p=gpd.GeoSeries(df1.geometry)
p2=gpd.GeoSeries(df2.geometry)
p3=gpd.GeoSeries(df3.geometry)
p4=gpd.GeoSeries(df4.geometry)

p.plot(alpha=0.4,color="blue")
p2.plot(alpha=0.4,color="red")
plt.scatter(Datos["Longitude"],Datos["Latitude"])


ax1=p.plot(alpha=0.3)
p2.plot(ax=ax1,alpha=0.3,color="red")
p3.plot(ax=ax1,alpha=0.3,color="fuchsia")
p4.plot(ax=ax1,alpha=0.3,color="black")
plt.scatter(Datos["Longitude"],Datos["Latitude"],s=np.exp(5*(M-np.min(M))/(np.max(M)-np.min(M))+1),alpha=0.5)
plt.xlabel("Longitud")
plt.ylabel("Latitud")
plt.title("Sismos")
plt.scatter(Datos["Longitude"],Datos["Latitude"],s=0.5)





#################


M=Datos["Magnitude"]



numdays=100





flag=0
flagAux=0
for i in range(len(Datos)):
    fig, (ax1, ax2) = plt.subplots(nrows=2,height_ratios=[10, 1])
    if GNSSDate[GNSSDate<Fechas[i]][-1] and (yder[GNSSDate<Fechas[i]][-1])<0:
        if flag!=flagAux:
            flag+=1
            flagAux=flag
            
        if flag==1:
            ax1=p.plot(ax=ax1,alpha=0.3)#p.plot(alpha=0.3)
            ax1.set_xlabel("Longitude")    
            ax1.set_ylabel("Latitude")    
            ax1.set_xlim(np.min(Datos["Longitude"]),np.max(Datos["Longitude"]))    
            ax1.set_ylim(np.min(Datos["Latitude"]),np.max(Datos["Latitude"]))    
            ax1.set_aspect('equal')        


        elif flag==2:
            ax1=p2.plot(ax=ax1,alpha=0.3)
            ax1.set_xlabel("Longitude")    
            ax1.set_ylabel("Latitude")    
            ax1.set_xlim(np.min(Datos["Longitude"]),np.max(Datos["Longitude"]))    
            ax1.set_ylim(np.min(Datos["Latitude"]),np.max(Datos["Latitude"]))    
            ax1.set_aspect('equal')        


        elif flag==3:
            ax1=p3.plot(ax=ax1,alpha=0.3)
            ax1.set_xlabel("Longitude")    
            ax1.set_ylabel("Latitude")    
            ax1.set_xlim(np.min(Datos["Longitude"]),np.max(Datos["Longitude"]))    
            ax1.set_ylim(np.min(Datos["Latitude"]),np.max(Datos["Latitude"]))    
            ax1.set_aspect('equal')        

        

        elif flag==4:
            1
            #ax1=p4.plot(ax=ax1,alpha=0.3)
        #plt.gca().set_aspect('equal')        

    
    else :
        flagAux=-1

    


    Select=(Fechas<Fechas[i])

    for j in range(len(Fechas[Select])):
        Select[j]=(Fechas[j]-Fechas[i]).total_seconds()>-sectoday*numdays

    Long=Datos["Longitude"][Select]
    Lat=Datos["Latitude"][Select]
    Size=5+np.exp(5*(M[Select]-np.min(M))/(np.max(M)-np.min(M))+1)
    
    
    ax1.scatter(Long,Lat,Size, alpha=0.5)
    ax1.set_xlim(np.min(Datos["Longitude"]),np.max(Datos["Longitude"]))    
    ax1.set_ylim(np.min(Datos["Latitude"]),np.max(Datos["Latitude"]))   
    ax1.set_xlabel("Longitude")    
    ax1.set_ylabel("Latitude")    
    ax1.set_aspect('equal')       
    ax1.plot(Trench["x"],Trench["y"],color="hotpink")
    #plt.gca().set_aspect('equal')
    ###############Branch structure
    PIJReduc=PIJ[Select,:][:,Select]
    for r in range(len(PIJReduc)):
        for s in range(len(PIJReduc)):
            if PIJReduc[s,r]>=0.2:
                ax1.annotate("", xy=(Long[Long.index[s]],Lat[Long.index[s]]), xytext=(Long[Long.index[r]],Lat[Long.index[r]]),arrowprops=dict(arrowstyle="->,head_width=0.1",alpha=PIJReduc[s,r]))
    
        
    shape=shp.Reader("./Mexico.shp") #http://geoportal.conabio.gob.mx/metadatos/doc/html/dest_2010gw.html
    
    
    for shape in shape.shapeRecords():
        x = [i[0] for i in shape.shape.points[:]]
        y = [i[1] for i in shape.shape.points[:]]
        ax1.plot(x,y,color="black")

    ################
    ax2.set_ylim(-1,1)
    ax2.plot(Fechas,np.zeros(len(Fechas)))
    ax2.set_xlabel("Date")    
    ax2.axvline(Fechas[i])
    ax2.yaxis.set_visible(False)
    

    fig.savefig("/Users/isaias/Desktop/SSE_GNSS/"+str(i)+".jpg", dpi=200)
    plt.show()    

        
import imageio
images = []
for i in range(len(Datos)):
    images.append(imageio.imread("/Users/isaias/Desktop/SSE_GNSS/"+str(i)+".jpg"))

imageio.mimsave('/Users/isaias/Desktop/movie.gif', images)


###############################
ind=np.arange(183,193)
# ind=Datos.index
LongGraf=Datos.loc[ind]["Longitude"]
LatGraf=Datos.loc[ind]["Latitude"]
M2=Datos.loc[ind]["Magnitude"]
Size=np.exp(5*(M2-np.min(M2))/(np.max(M2)-np.min(M2))+1)

# shape=shp.Reader("./Mexico.shp") #http://geoportal.conabio.gob.mx/metadatos/doc/html/dest_2010gw.html

# for shape in shape.shapeRecords():
#     x = [i[0] for i in shape.shape.points[:]]
#     y = [i[1] for i in shape.shape.points[:]]
#     plt.plot(x,y,color="black",alpha=0.3)
plt.xlim(np.min(LongGraf-0.1),np.max(LongGraf+0.1))
plt.ylim(np.min(LatGraf-0.1),np.max(LatGraf+0.1))
# plt.plot(Trench["x"],Trench["y"],color="hotpink")
# plt.scatter(Reference[:,0],Reference[:,1],s=np.exp(5*(MRef-np.min(M))/(np.max(M)-np.min(M))+1),alpha=0.5,color="Green")
# for i in range(len(Reference)):
#     plt.text(Reference[i,0]-0.2,Reference[i,1]+0.2, str(i+4.5))

plt.scatter(LongGraf,LatGraf,s=Size,alpha=0.5)
for j in ind:
    for s in ind:
        if PIJ[s,j]>0.2:
            plt.annotate("", xy=Datos[["Longitude","Latitude"]].loc[s], xytext=Datos[["Longitude","Latitude"]].loc[j],arrowprops=dict(arrowstyle="->,head_width=0.1",alpha=PIJ[s,j]))
            print(j,s,PIJ[s,j])

plt.xlabel("Longitude")
plt.ylabel("Latitude")
PIJ[192,183]
PIJ[s,j]




################################
ax1=p.plot(alpha=0.3)
plt.xlim(-102.5,-96)
plt.ylim(15,19.5)
plt.plot(Trench["x"],Trench["y"],color="hotpink")
plt.ylabel("Latitude")
plt.xlabel("Longitude")

shape=shp.Reader("./Mexico.shp") #http://geoportal.conabio.gob.mx/metadatos/doc/html/dest_2010gw.html


plt.xlim(DomX[0],DomX[-1])
plt.ylim(DomY[0],DomY[-1])
for shape in shape.shapeRecords():
    x = [i[0] for i in shape.shape.points[:]]
    y = [i[1] for i in shape.shape.points[:]]
    plt.plot(x,y,color="black")







plt.xlim(-102.5,-96)
plt.ylim(15,19.5)
ax1=p2.plot(alpha=0.3)
plt.plot(Trench["x"],Trench["y"],color="hotpink")
plt.xlim(-102.5,-96)
plt.ylim(15,19.5)
plt.ylabel("Latitude")
plt.xlabel("Longitude")

shape=shp.Reader("./Mexico.shp") #http://geoportal.conabio.gob.mx/metadatos/doc/html/dest_2010gw.html


plt.xlim(DomX[0],DomX[-1])
plt.ylim(DomY[0],DomY[-1])
for shape in shape.shapeRecords():
    x = [i[0] for i in shape.shape.points[:]]
    y = [i[1] for i in shape.shape.points[:]]
    plt.plot(x,y,color="black")



plt.xlim(-102.5,-96)
plt.ylim(15,19.5)
ax1=p3.plot(alpha=0.3)
plt.xlim(-102.5,-96)
plt.ylim(15,19.5)
plt.plot(Trench["x"],Trench["y"],color="hotpink")
plt.ylabel("Latitude")
plt.xlabel("Longitude")

shape=shp.Reader("./Mexico.shp") #http://geoportal.conabio.gob.mx/metadatos/doc/html/dest_2010gw.html
plt.xlim(DomX[0],DomX[-1])
plt.ylim(DomY[0],DomY[-1])
for shape in shape.shapeRecords():
    x = [i[0] for i in shape.shape.points[:]]
    y = [i[1] for i in shape.shape.points[:]]
    plt.plot(x,y,color="black")




################################

MAX=np.max(np.hstack((mu0,mu1,mu2,mu3))) #Falta mu4    
MIN=np.min(np.hstack((mu0,mu1,mu2,mu3))) #Falta mu4    
# plt.imshow(mu0[:,::-1], extent=[DomX[0],DomX[1],DomY[0],DomY[1]],aspect=DeltaX/DeltaY)

plt.title('$\mu^0$')


              

plt.imshow(mu0[::-1,:], extent=[DomX[0],DomX[-1],DomY[0],DomY[-1]],aspect=DeltaX/DeltaY,vmin=np.min(MIN,0),vmax=MAX)
plt.colorbar()
plt.ylabel("Latitude")
plt.plot(Trench["x"],Trench["y"],color="hotpink")
plt.xlim(DomX[0],DomX[-1])
plt.ylim(DomY[0],DomY[-1])
#plt.scatter(Datos["Longitude"],Datos["Latitude"],s=np.exp(5*(M-np.min(M))/(np.max(M)-np.min(M))+1),alpha=0.1)




ax1=p.plot(facecolor="none", edgecolor='white')
ax1.set_title('$\mu^1$')
ax1.imshow(mu1[::-1,:], extent=[DomX[0],DomX[-1],DomY[0],DomY[-1]],aspect=DeltaX/DeltaY,vmin=0,vmax=MAX)    
plt.plot(Trench["x"],Trench["y"],color="hotpink")
plt.xlim(DomX[0],DomX[-1])
plt.ylim(DomY[0],DomY[-1])



ax2=p2.plot(facecolor="none", edgecolor='white')
ax2.set_title('$\mu^2$')
ax2.imshow(mu2[::-1,:], extent=[DomX[0],DomX[-1],DomY[0],DomY[-1]],aspect=DeltaX/DeltaY,vmin=0,vmax=MAX)
plt.plot(Trench["x"],Trench["y"],color="hotpink")
plt.xlim(DomX[0],DomX[-1])
plt.ylim(DomY[0],DomY[-1])
plt.ylabel("Latitude")
plt.xlabel("Longitude")



ax3=p3.plot(facecolor="none", edgecolor='white')
ax3.set_title('$\mu^3$')
ax3.imshow(mu3[::-1,:], extent=[DomX[0],DomX[-1],DomY[0],DomY[-1]],aspect=DeltaX/DeltaY,vmin=0,vmax=MAX)
plt.plot(Trench["x"],Trench["y"],color="hotpink")
plt.xlim(DomX[0],DomX[-1])
plt.ylim(DomY[0],DomY[-1])
plt.xlabel("Longitude")






#################################
fig3, axs = plt.subplots(nrows=1, ncols=1, figsize=(14, 8))
axs.scatter(np.arange(len(PIJ)),np.apply_along_axis(sum, 0,PIJ))
#axs.set_title(thetaopt)
plt.xlabel("Index")
plt.ylabel("Expected Aftershocks")
################################

plt.hist(np.apply_along_axis(sum, 1,PII))
plt.xlabel('$\sum_{s=0}^3 {\hat{p}^s_{ii}}$')    
plt.ylabel('Frequency')    


    

###############################
















    



    
    
    
    

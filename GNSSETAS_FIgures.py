#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 11:47:05 2024

@author: isaias
"""

#https://www.gsi.go.jp/kankyochiri/gm_japan_e.html

#Run GNSSETAS_EM2 before
import shapefile as shp
from mpmath import mp
import numpy as np
import pandas as pd
import scipy as sp
import geopy.distance
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import time
import shapefile as shp
from multiprocessing import Pool, cpu_count, get_context, Process

import pwlf
#################################


Earthquakes2=pd.read_csv('./Earthquakes2.csv')
xlim,Xlim=lonmin,lonmax
ylim,Ylim=latmin,latmax
Earthquakes2=Earthquakes2[~Earthquakes2["Strike 2"].isnull()]
Earthquakes2=Earthquakes2.loc[ (Earthquakes2["Longitude"]<Xlim) & (Earthquakes2["Longitude"]>xlim) &  (Earthquakes2["Latitude"]<Ylim)&  (Earthquakes2["Latitude"]>ylim) ]
#Earthquakes=Earthquakes.loc[Earthquakes["Depth"]!="Sh"]
#Earthquakes["Depth"]=Earthquakes["Depth"].astype(float) 
Earthquakes2=Earthquakes2.loc[Earthquakes2["Depth"]<=100]
Earthquakes2=Earthquakes2.loc[Earthquakes2["Mw*"]>=M0]
# Earthquakes2=Earthquakes2.loc[Earthquakes2["Year"]>=2000]
Earthquakes2=Earthquakes2.loc[Earthquakes2["Year"]<=anio]
Earthquakes2=Earthquakes2.loc[Earthquakes2["Year"]>=2001]
# Earthquakes=Earthquakes.loc[(Earthquakes["Rake 1"]>45) & (Earthquakes["Rake 1"]<134)]
Earthquakes2=Earthquakes2[Earthquakes2.apply(inside, axis=1)]
Earthquakes2=Earthquakes2.reset_index(drop=True)




sectoday=86400
for i in range(len(Earthquakes2)):
  if i==0:
    Fechas2=datetime(Earthquakes2["Year"][i], Earthquakes2["Month"][i], Earthquakes2["Day"][i], int(Earthquakes2["Hour"][i]), int(Earthquakes2["Minute"][i]), int(Earthquakes2["Second"][i]), int(str(Earthquakes2["Second"][i])[-1])*100000)
  else :
    Fechas2=np.hstack((Fechas2,datetime(Earthquakes2["Year"][i], Earthquakes2["Month"][i], Earthquakes2["Day"][i], int(Earthquakes2["Hour"][i]), int(Earthquakes2["Minute"][i]), int(Earthquakes2["Second"][i]), int(str(Earthquakes2["Second"][i])[-1])*100000)))


Diferencias2=np.zeros(len(Fechas2))
for i in range(len(Diferencias2)):
    Diferencias2[i]=(Fechas2[i]-mindate).total_seconds()

x = np.array(Diferencias2)
y = np.arange(len(Fechas2))

my_pwlf = pwlf.PiecewiseLinFit(x, y)
breaks = my_pwlf.fit(5)
print(breaks)

fechagrafica=np.copy(Fechas2)
for i in range(len(Diferencias2)):
    fechagrafica[i]=mindate+timedelta(seconds=Diferencias2[i])

cortes=np.copy(Fechas2[:len(breaks)])
for i in range(len(cortes)):
    cortes[i]=mindate+timedelta(seconds=breaks[i])


y_hat = my_pwlf.predict(x)

plt.figure()
plt.plot(fechagrafica, y, 'o')
plt.plot(fechagrafica, y_hat, '-')
plt.show()


#################################

# plt.plot(Fechas,np.arange(len(Fechas)),label="raw")

# plt.plot(Fechas2,np.arange(len(Fechas2)),label="Sawires")
# plt.plot(fechagrafica, y_hat, '-',linewidth=1,alpha=0.5)
# for i in range(len(cortes)):
#     plt.axvline(cortes[i],linewidth=1,color="purple",alpha=0.5)
# plt.show()

for i in range(len(Fechas[Datos["Magnitude"]>7])):
    plt.axvline(Fechas[Datos["Magnitude"]>7][i],linewidth=0.5,color="cyan",alpha=0.5)

for i in range(len(DominioInt)):
    plt.axvline(GNSSDate[DominioInt][i],color="firebrick")


Extension='./figures/'
PII=np.loadtxt(Extension+'PII.csv', delimiter=',')
background=np.apply_along_axis(np.sum,1,PII)>0.99


plt.plot(Fechas[background &  (Fechas>datetime(2001, 1, 1))],np.arange(len(Fechas[background &  (Fechas>datetime(2001, 1, 1))])),label="Full")
for i in range(len(Fechas[Datos["Magnitude"]>7])):
    plt.axvline(Fechas[Datos["Magnitude"]>7][i],linewidth=0.5,color="black",alpha=0.5)
plt.show()

# 
PIINull=np.loadtxt(Extension+'PIINull.csv', delimiter=',')
len(PIINull)
background2=np.apply_along_axis(np.sum,1,PIINull)>0.99
np.sum(background2)
# plt.plot(Fechas[background2 &  (Fechas>datetime(2001, 1, 1))],np.arange(len(Fechas[background2 &  (Fechas>datetime(2001, 1, 1))])),label="Reduced")
for i in range(len(Fechas[Datos["Magnitude"]>7])):
    plt.axvline(Fechas[Datos["Magnitude"]>7][i],linewidth=0.5,color="black",alpha=0.5)
plt.xlabel("Date")
plt.ylabel("N(t)")
plt.legend()
plt.show()


################
################




import geopandas as gpd
# Trench=pd.read_csv('./trench_cocos.csv', delimiter=',')
# M=Datos["Magnitude"]




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



file_path = './Cocos_Plate_Z120.dat'
df = pd.read_csv(file_path, delimiter='\t')
df=df.loc[df[df.columns[0]]>DomX[0]] 
df=df.loc[df[df.columns[0]]<DomX[-1]]
df=df.loc[df[df.columns[1]]>DomY[0]]
df=df.loc[df[df.columns[1]]<DomY[-1]] 
np.min(df[df.columns[1]])
Xslab, Yslab = np.meshgrid(np.arange(lonmin,lonmax,0.01), np.arange(latmin,latmax,0.01))
Zslab = griddata((df[df.columns[0]], df[df.columns[1]]), df[df.columns[2]], (Xslab, Yslab), method='cubic')
np.max(df[df.columns[2]])



shape=shp.Reader("./Mexico.shp") #http://geoportal.conabio.gob.mx/metadatos/doc/html/dest_2010gw.html

Trench=pd.read_csv('./trench_cocos.csv', delimiter=',')
Reference=np.array([[-104.5,16],[-104,16],[-103.5,16],[-103,16]])
MRef=np.array([4.5,5.5,6.5,7.5])


fig, ax = plt.subplots()

fig.gca().set_facecolor('lightblue')
ax.set_xlim(DomX[0],DomX[-1])
ax.set_ylim(DomY[0],DomY[-1])
for shape in shape.shapeRecords():
    x = [i[0] for i in shape.shape.points[:]]
    y = [i[1] for i in shape.shape.points[:]]
    # plt.plot(x,y,color="gray")
    plt.fill(x, y, color="yellowgreen", edgecolor="lightgray")


contours=plt.contour(Xslab, Yslab, -Zslab,levels=4,vmin=0, colors='darkviolet')
plt.clabel(contours, inline=True, fontsize=12)
# plt.colorbar()   

gdf = gpd.read_file("FFM2003.geojson")  # Replace with your file path
# gdf.plot(ax=ax,column='slip',edgecolor='black', cmap='hot',figsize=(10, 10))

centroids = gdf.geometry.centroid
x = centroids.x.values
y = centroids.y.values
z = gdf['slip'].values  # Replace with your actual data column
# Create a grid for interpolation
xi = np.linspace(x.min(), x.max(), 200)
yi = np.linspace(y.min(), y.max(), 200)
xi, yi = np.meshgrid(xi, yi)
zi = griddata((x, y), z, (xi, yi), method='cubic')
contours = ax.contour(xi, yi, zi, levels=np.arange(1,6.5,1), cmap='hot',alpha=0.5)
# ax.clabel(contours, inline=True, fontsize=8)
# plt.show()


# 
gdf2 = gpd.read_file("FFM2012.geojson")  # Replace with your file path
centroids = gdf2.geometry.centroid
x = centroids.x.values
y = centroids.y.values
z = gdf2['slip'].values  # Replace with your actual data column
# Create a grid for interpolation
xi = np.linspace(x.min(), x.max(), 200)
yi = np.linspace(y.min(), y.max(), 200)
xi, yi = np.meshgrid(xi, yi)
zi = griddata((x, y), z, (xi, yi), method='cubic')
contours = ax.contour(xi, yi, zi, levels=np.arange(1,6.5,1), cmap='hot',alpha=0.5)
# 
gdf3 = gpd.read_file("FFM2014.geojson")  # Replace with your file path
centroids = gdf3.geometry.centroid
x = centroids.x.values
y = centroids.y.values
z = gdf3['slip'].values  # Replace with your actual data column
# Create a grid for interpolation
xi = np.linspace(x.min(), x.max(), 200)
yi = np.linspace(y.min(), y.max(), 200)
xi, yi = np.meshgrid(xi, yi)
zi = griddata((x, y), z, (xi, yi), method='cubic')
contours = ax.contour(xi, yi, zi, levels=np.arange(1,6.5,1), cmap='hot',alpha=0.5)

##################################################

p=gpd.GeoSeries(df1.geometry[4])
p2=gpd.GeoSeries(df2.geometry[3])
p3=gpd.GeoSeries(df3.geometry[3])
p4=gpd.GeoSeries(df4.geometry)

p.plot(ax=ax,alpha=0.8,facecolor='none',edgecolor="blue")

p2.plot(ax=ax,alpha=0.8,edgecolor="orange",facecolor='none')
p3.plot(ax=ax,alpha=0.8,edgecolor="fuchsia",facecolor='none')
p4.plot(ax=ax,alpha=0.8,edgecolor="black",facecolor='none')
##################################################




for i in range(len(DomX)):
    plt.axvline(DomX[i],color="firebrick",alpha=0.5,linewidth=0.5)
    plt.axhline(DomY[i],color="firebrick",alpha=0.5,linewidth=0.5)
plt.plot(Trench["x"],Trench["y"],color="hotpink")
plt.scatter(Datos["Longitude"],Datos["Latitude"],s=np.exp(5*(M-np.min(M))/(np.max(M)-np.min(M))+1),alpha=0.5,color="dodgerblue")
# plt.scatter(Earthquakes2["Longitude"],Earthquakes2["Latitude"],s=np.exp(5*(Earthquakes2["Mw*"]-np.min(Earthquakes2["Mw*"]))/(np.max(Earthquakes2["Mw*"])-np.min(Earthquakes2["Mw*"]))+1),alpha=0.5,color="dodgerblue")


plt.scatter(Reference[:,0],Reference[:,1],s=np.exp(5*(MRef-np.min(M))/(np.max(M)-np.min(M))+1),alpha=0.5,color="Green")
for i in range(len(Reference)):
    plt.text(Reference[i,0]-0.2,Reference[i,1]+0.2, str(i+4.5))

plt.ylabel("Latitude")
plt.xlabel("Longitude")


plt.scatter(-100.2672,17.0485, color="gold",marker="x")
plt.show()


################################################################################

# gdf.plot(column="slip")
# plt.tight_layout()
# plt.show()

# import geopy.distance
# coords_1 = (-104.125, 18.8)
# coords_2 = (-104.125, 19.25)
# coords_1 = (18.8,-104.125)
# coords_2 = (19.25,-104.125)

# print(geopy.distance.geodesic(coords_1, coords_2).km)

################################################################################



    
################ Hyperparameter
Trench=pd.read_csv('./trench_cocos.csv', delimiter=',')
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

# # Datos=pd.read_csv('/Users/isaias/Desktop/hypo_MASUB.csv',delimiter=",")
# Datos=pd.read_csv('./Earthquakes.csv',delimiter=",")
# Datos=Datos.loc[ Datos["Year"]<anio]
# M0=3.5
# Datos=Datos.loc[ Datos["Magnitude"]>=M0]
# Datos=Datos.reset_index()



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
# kde =KernelDensity(kernel="gaussian",bandwidth=10_000_000).fit(np.reshape(Diferencias[Selected],(len(Diferencias[Selected]),1)))    
# plt.plot(GNSSDif,GNSSData[:,3])        
# # plt.show()
# plt.plot(Diferencias[Selected],10**(7)*np.exp(kde.score_samples(np.reshape(Diferencias[Selected],(len(Diferencias[Selected]),1))))) 
# plt.show()    

# fig = px.line(x=GNSSDif/sectoday, y=GNSSData[:,3] )
# #fig.add_trace(go.Scatter(x=Diferencias, y=Datos["Magnitude"], mode='markers'))
# fig.add_trace(go.Scatter(x=Diferencias[Selected]/sectoday,y=10**(7)*np.exp(kde.score_samples(np.reshape(Diferencias[Selected],(len(Diferencias[Selected]),1)))), marker_size=5*np.exp(Datos[Selected]["Magnitude"]-4), mode='lines+markers'))
# fig.show()
        

# plt.hist(Datos[Selected]["Magnitude"])
# plt.show() 
    
    










import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

x=GNSSDif/sectoday
y=GNSSData[:,3] 
# y=sp.signal.detrend(y) #Detrend
plt.plot(GNSSDate,y)
plt.ylabel("Desplazamiento [cm]")
tck = interpolate.splrep(x, y, s=0.05)
plt.show()

xnew = x
ynew = interpolate.splev(xnew, tck, der=0)
plt.plot(x,y)
plt.plot(xnew,ynew)
yder = interpolate.splev(xnew, tck, der=1)   # or BSpline(*tck)(xnew, 1)
plt.plot(x,yder)
plt.show()


plt.scatter(xnew,ynew)
plt.scatter(xnew[yder<0],(ynew[yder<0]))
plt.scatter(xnew[yder<0],100*np.abs(yder[yder<0]))
plt.show()


cumsum=np.cumsum(yder<0)
plt.plot(cumsum)
DominioInt=np.zeros(0,dtype=int)
plt.show()


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
plt.show()



    
#############################    
    
    

import geopandas as gpd
Trench=pd.read_csv('./trench_cocos.csv', delimiter=',')
M=Datos["Magnitude"]


gdf = gpd.read_file("FFM2003.geojson")  # Replace with your file path
gdf.slip
gdf.plot(column='slip',edgecolor='black', figsize=(10, 10))
plt.show()



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
plt.show()

ax1=p.plot(alpha=0.3)
p2.plot(ax=ax1,alpha=0.3,color="red")
p3.plot(ax=ax1,alpha=0.3,color="fuchsia")
p4.plot(ax=ax1,alpha=0.3,color="black")
plt.scatter(Datos["Longitude"],Datos["Latitude"],s=np.exp(5*(M-np.min(M))/(np.max(M)-np.min(M))+1),alpha=0.5)
plt.xlabel("Longitud")
plt.ylabel("Latitud")
plt.title("Sismos")
plt.scatter(Datos["Longitude"],Datos["Latitude"],s=0.5)
plt.show()

    
    

#################
import shapefile as shp

M=Datos["Magnitude"]



numdays=100



i=1473

flag=0
flagAux=0
for i in np.arange(0,len(Datos)):
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
            # 1
            ax1=p4.plot(ax=ax1,alpha=0.3)
            ax1.set_xlabel("Longitude")    
            ax1.set_ylabel("Latitude")    
            ax1.set_xlim(np.min(Datos["Longitude"]),np.max(Datos["Longitude"]))    
            ax1.set_ylim(np.min(Datos["Latitude"]),np.max(Datos["Latitude"]))    
            ax1.set_aspect('equal')        


    
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
        if PIJ[s,j]>0.1:
            plt.annotate("", xy=Datos[["Longitude","Latitude"]].loc[s], xytext=Datos[["Longitude","Latitude"]].loc[j],arrowprops=dict(arrowstyle="->,head_width=0.1",alpha=PIJ[s,j]))
            print(j,s,PIJ[s,j])

plt.xlabel("Longitude")
plt.ylabel("Latitude")
PIJ[192,183]
PIJ[s,j]


for s in np.arange(182,193):
    plt.annotate(xy=Datos[["Longitude","Latitude"]].loc[s], text=str(Datos.index[s]))#Datos["index"].loc[s]))
    

plt.show()
################################
Trench=pd.read_csv('./trench_cocos.csv', delimiter=',')
import shapefile as shp

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

plt.show()


################################




shape=shp.Reader("./Mexico.shp") #http://geoportal.conabio.gob.mx/metadatos/doc/html/dest_2010gw.html





import geopandas as gpd
Trench=pd.read_csv('./trench_cocos.csv', delimiter=',')
M=Datos["Magnitude"]


gdf = gpd.read_file("FFM2003.geojson")  # Replace with your file path
gdf.slip
gdf.plot(column='slip',edgecolor='black', figsize=(10, 10))
plt.show()



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


MAX=np.max(np.hstack((mu0,mu1,mu2,mu3,mu4))) #Falta mu4    
MIN=np.min(np.hstack((mu0,mu1,mu2,mu3,mu4))) #Falta mu4    
# plt.imshow(mu0[:,::-1], extent=[DomX[0],DomX[1],DomY[0],DomY[1]],aspect=DeltaX/DeltaY)





fig, ax = plt.subplots()
# ax.set_title('$\mu^0$')
ax.set_xlim(DomX[0],DomX[-1])
ax.set_ylim(DomY[0],DomY[-1])
for shape2 in shape.shapeRecords():
    x = [i[0] for i in shape2.shape.points[:]]
    y = [i[1] for i in shape2.shape.points[:]]
    # plt.plot(x,y,color="gray")
    plt.fill(x, y, color="none", edgecolor="gray")

# plt.show()

# im=ax.imshow(mu0[::-1,:], extent=[DomX[0],DomX[-1],DomY[0],DomY[-1]],aspect=0.5*DeltaX/DeltaY,vmin=np.min(MIN,0),vmax=MAX)
im=ax.imshow(mu0[::-1,:], extent=[DomX[0],DomX[-1],DomY[0],DomY[-1]],aspect=0.5*DeltaX/DeltaY)#,vmin=np.min(0),vmax=0.035)
cbar = fig.colorbar(im, ax=ax,fraction=0.022, pad=0.05)
# cbar.ax.set_title(r'$\frac{\text{Events}}{\text{deg}^2 \text{day}}$', fontsize=12, pad=10)
cbar.set_label(r'$\frac{\text{Events}}{\text{deg}^2 \text{day}}$', fontsize=12)#, pad=10)
ax.set_ylabel("Latitude")
ax.set_xlabel("Longitude")
ax.plot(Trench["x"],Trench["y"],color="hotpink")
ax.set_xlim(DomX[0],DomX[-1])
ax.set_ylim(DomY[0],DomY[-1])
#plt.scatter(Datos["Longitude"],Datos["Latitude"],s=np.exp(5*(M-np.min(M))/(np.max(M)-np.min(M))+1),alpha=0.1)
plt.show()




fig, ax = plt.subplots()
# ax.set_title('$\mu^1$')
ax.set_xlim(DomX[0],DomX[-1])
ax.set_ylim(DomY[0],DomY[-1])
for shape2 in shape.shapeRecords():
    x = [i[0] for i in shape2.shape.points[:]]
    y = [i[1] for i in shape2.shape.points[:]]
    # plt.plot(x,y,color="gray")
    plt.fill(x, y, color="none", edgecolor="gray")

im=ax.imshow(mu1[::-1,:], extent=[DomX[0],DomX[-1],DomY[0],DomY[-1]],aspect=0.5*DeltaX/DeltaY,vmin=np.min(MIN,0),vmax=MAX)
# im=ax.imshow(mu1[::-1,:], extent=[DomX[0],DomX[-1],DomY[0],DomY[-1]],aspect=0.5*DeltaX/DeltaY)
cbar =fig.colorbar(im, ax=ax,fraction=0.022, pad=0.05)
cbar.set_label(r'$\frac{\text{Events}}{\text{deg}^2 \text{day}}$', fontsize=12)#, pad=10)
# im=ax.imshow(mu1[::-1,:], extent=[DomX[0],DomX[-1],DomY[0],DomY[-1]],aspect=0.5*DeltaX/DeltaY,vmin=np.min(MIN,0),vmax=MAX)
p.plot(ax=ax,facecolor="none", edgecolor='white')
ax.set_ylabel("Latitude")
ax.set_xlabel("Longitude")
ax.plot(Trench["x"],Trench["y"],color="hotpink")
ax.set_xlim(DomX[0],DomX[-1])
ax.set_ylim(DomY[0],DomY[-1])
#plt.scatter(Datos["Longitude"],Datos["Latitude"],s=np.exp(5*(M-np.min(M))/(np.max(M)-np.min(M))+1),alpha=0.1)
plt.show()



fig, ax = plt.subplots()
# ax.set_title('$\mu^2$')
ax.set_xlim(DomX[0],DomX[-1])
ax.set_ylim(DomY[0],DomY[-1])
for shape2 in shape.shapeRecords():
    x = [i[0] for i in shape2.shape.points[:]]
    y = [i[1] for i in shape2.shape.points[:]]
    # plt.plot(x,y,color="gray")
    plt.fill(x, y, color="none", edgecolor="gray")

# plt.imshow(mu0[::-1,:], extent=[DomX[0],DomX[-1],DomY[0],DomY[-1]],aspect=DeltaX/DeltaY,vmin=np.min(MIN,0),vmax=MAX)
im=ax.imshow(mu2[::-1,:], extent=[DomX[0],DomX[-1],DomY[0],DomY[-1]],aspect=0.5*DeltaX/DeltaY,vmin=np.min(MIN,0),vmax=MAX)
p2.plot(ax=ax,facecolor="none", edgecolor='white')

# im=ax.imshow(mu2[::-1,:], extent=[DomX[0],DomX[-1],DomY[0],DomY[-1]],aspect=0.5*DeltaX/DeltaY)
cbar = fig.colorbar(im, ax=ax,fraction=0.022, pad=0.05)
cbar.set_label(r'$\frac{\text{Events}}{\text{deg}^2 \text{day}}$', fontsize=12)#, pad=10)p2.plot(ax=ax,facecolor="none", edgecolor='white')
# fig.colorbar(im, ax=ax,fraction=0.03, pad=0.01)
ax.set_ylabel("Latitude")
ax.set_xlabel("Longitude")
ax.plot(Trench["x"],Trench["y"],color="hotpink")
ax.set_xlim(DomX[0],DomX[-1])
ax.set_ylim(DomY[0],DomY[-1])
#plt.scatter(Datos["Longitude"],Datos["Latitude"],s=np.exp(5*(M-np.min(M))/(np.max(M)-np.min(M))+1),alpha=0.1)
gdf = gpd.read_file("FFM2003.geojson")  # Replace with your file path
centroids = gdf.geometry.centroid
x = centroids.x.values
y = centroids.y.values
z = gdf['slip'].values  # Replace with your actual data column
# Create a grid for interpolation
xi = np.linspace(x.min(), x.max(), 200)
yi = np.linspace(y.min(), y.max(), 200)
xi, yi = np.meshgrid(xi, yi)
zi = griddata((x, y), z, (xi, yi), method='cubic')
contours = ax.contour(xi, yi, zi, levels=np.arange(1,6.5,1), cmap='hot',alpha=0.5)
# ax.clabel(contours, inline=True, fontsize=8)
# plt.show()
plt.show()




fig, ax = plt.subplots()
# ax.set_title('$\mu^3$')
ax.set_xlim(DomX[0],DomX[-1])
ax.set_ylim(DomY[0],DomY[-1])
for shape2 in shape.shapeRecords():
    x = [i[0] for i in shape2.shape.points[:]]
    y = [i[1] for i in shape2.shape.points[:]]
    # plt.plot(x,y,color="gray")
    plt.fill(x, y, color="none", edgecolor="gray")
im=ax.imshow(mu3[::-1,:], extent=[DomX[0],DomX[-1],DomY[0],DomY[-1]],aspect=0.5*DeltaX/DeltaY,vmin=np.min(MIN,0),vmax=MAX)
# im=ax.imshow(mu3[::-1,:], extent=[DomX[0],DomX[-1],DomY[0],DomY[-1]],aspect=0.5*DeltaX/DeltaY)    
cbar = fig.colorbar(im, ax=ax,fraction=0.022, pad=0.05)
cbar.set_label(r'$\frac{\text{Events}}{\text{deg}^2 \text{day}}$', fontsize=12)#, pad=10)
# cbar.ax.set_title(r'$\frac{\text{Events}}{\text{deg}^2 \text{day}}$', fontsize=12, pad=10)
# plt.imshow(mu0[::-1,:], extent=[DomX[0],DomX[-1],DomY[0],DomY[-1]],aspect=DeltaX/DeltaY,vmin=np.min(MIN,0),vmax=MAX)
p3.plot(ax=ax,facecolor="none", edgecolor='white')
# fig.colorbar(im, ax=ax,fraction=0.03, pad=0.01)
ax.set_ylabel("Latitude")
ax.set_xlabel("Longitude")
ax.plot(Trench["x"],Trench["y"],color="hotpink")
ax.set_xlim(DomX[0],DomX[-1])
ax.set_ylim(DomY[0],DomY[-1])
#plt.scatter(Datos["Longitude"],Datos["Latitude"],s=np.exp(5*(M-np.min(M))/(np.max(M)-np.min(M))+1),alpha=0.1)
gdf = gpd.read_file("FFM2003.geojson")  # Replace with your file path
centroids = gdf.geometry.centroid
x = centroids.x.values
y = centroids.y.values
z = gdf['slip'].values  # Replace with your actual data column
# Create a grid for interpolation
xi = np.linspace(x.min(), x.max(), 200)
yi = np.linspace(y.min(), y.max(), 200)
xi, yi = np.meshgrid(xi, yi)
zi = griddata((x, y), z, (xi, yi), method='cubic')
contours = ax.contour(xi, yi, zi, levels=np.arange(1,6.5,1), cmap='hot',alpha=0.5)
# ax.clabel(contours, inline=True, fontsize=8)
# plt.show()
plt.show()






fig, ax = plt.subplots()
# ax.set_title('$\mu^4$')

ax.set_xlim(DomX[0],DomX[-1])
ax.set_ylim(DomY[0],DomY[-1])
for shape2 in shape.shapeRecords():
    x = [i[0] for i in shape2.shape.points[:]]
    y = [i[1] for i in shape2.shape.points[:]]
    # plt.plot(x,y,color="gray")
    plt.fill(x, y, color="none", edgecolor="gray")
im=ax.imshow(mu4[::-1,:], extent=[DomX[0],DomX[-1],DomY[0],DomY[-1]],aspect=0.5*DeltaX/DeltaY,vmin=np.min(MIN,0),vmax=MAX)
# im=ax.imshow(mu4[::-1,:], extent=[DomX[0],DomX[-1],DomY[0],DomY[-1]],aspect=0.5*DeltaX/DeltaY)
cbar = fig.colorbar(im, ax=ax,fraction=0.022, pad=0.05)
cbar.set_label(r'$\frac{\text{Events}}{\text{deg}^2 \text{day}}$', fontsize=12)#, pad=10)
# plt.imshow(mu0[::-1,:], extent=[DomX[0],DomX[-1],DomY[0],DomY[-1]],aspect=DeltaX/DeltaY,vmin=np.min(MIN,0),vmax=MAX)
p4.plot(ax=ax,facecolor="none", edgecolor='white')
# fig.colorbar(im, ax=ax,fraction=0.03, pad=0.01)
ax.set_ylabel("Latitude")
ax.set_xlabel("Longitude")
ax.plot(Trench["x"],Trench["y"],color="hotpink")
ax.set_xlim(DomX[0],DomX[-1])
ax.set_ylim(DomY[0],DomY[-1])
#plt.scatter(Datos["Longitude"],Datos["Latitude"],s=np.exp(5*(M-np.min(M))/(np.max(M)-np.min(M))+1),alpha=0.1)
gdf = gpd.read_file("FFM2003.geojson")  # Replace with your file path
centroids = gdf.geometry.centroid
x = centroids.x.values
y = centroids.y.values
z = gdf['slip'].values  # Replace with your actual data column
# Create a grid for interpolation
xi = np.linspace(x.min(), x.max(), 200)
yi = np.linspace(y.min(), y.max(), 200)
xi, yi = np.meshgrid(xi, yi)
zi = griddata((x, y), z, (xi, yi), method='cubic')
contours = ax.contour(xi, yi, zi, levels=np.arange(1,6.5,1), cmap='hot',alpha=0.5)
# ax.clabel(contours, inline=True, fontsize=8)
# plt.show()
# 
gdf2 = gpd.read_file("FFM2012.geojson")  # Replace with your file path
centroids = gdf2.geometry.centroid
x = centroids.x.values
y = centroids.y.values
z = gdf2['slip'].values  # Replace with your actual data column
# Create a grid for interpolation
xi = np.linspace(x.min(), x.max(), 200)
yi = np.linspace(y.min(), y.max(), 200)
xi, yi = np.meshgrid(xi, yi)
zi = griddata((x, y), z, (xi, yi), method='cubic')
contours = ax.contour(xi, yi, zi, levels=np.arange(1,6.5,1), cmap='hot',alpha=0.5)
# 
gdf3 = gpd.read_file("FFM2014.geojson")  # Replace with your file path
centroids = gdf3.geometry.centroid
x = centroids.x.values
y = centroids.y.values
z = gdf3['slip'].values  # Replace with your actual data column
# Create a grid for interpolation
xi = np.linspace(x.min(), x.max(), 200)
yi = np.linspace(y.min(), y.max(), 200)
xi, yi = np.meshgrid(xi, yi)
zi = griddata((x, y), z, (xi, yi), method='cubic')
contours = ax.contour(xi, yi, zi, levels=np.arange(1,6.5,1), cmap='hot',alpha=0.5)

plt.show()







#################################
fig3, axs = plt.subplots(nrows=1, ncols=1, figsize=(14, 8))
axs.scatter(np.arange(len(PIJ)),np.apply_along_axis(sum, 0,PIJ))
#axs.set_title(thetaopt)
plt.xlabel("Index")
plt.ylabel("Expected Aftershocks")
plt.show()
################################

plt.hist(np.apply_along_axis(sum, 1,PII))
plt.xlabel('$\sum_{s=0}^4 {\hat{p}^s_{ii}}$')    
plt.ylabel('Frequency')    

plt.show()
    

###############################



Trench=pd.read_csv('./trench_cocos.csv', delimiter=',')

fig, ax = plt.subplots()
# plt.imshow(mu0[::-1,:], extent=[DomX[0],DomX[-1],DomY[0],DomY[-1]],aspect=DeltaX/DeltaY,vmin=np.min(MIN,0),vmax=MAX)
im=ax.imshow(mu0[::-1,:], extent=[DomX[0],DomX[-1],DomY[0],DomY[-1]],aspect=0.5*DeltaX/DeltaY)
fig.colorbar(im, ax=ax,fraction=0.03, pad=0.01)
ax.set_ylabel("Latitude")
ax.plot(Trench["x"],Trench["y"],color="hotpink")
ax.set_xlim(DomX[0],DomX[-1])
ax.set_ylim(DomY[0],DomY[-1])
#plt.scatter(Datos["Longitude"],Datos["Latitude"],s=np.exp(5*(M-np.min(M))/(np.max(M)-np.min(M))+1),alpha=0.1)
plt.show()



ax1=p.plot(facecolor="none", edgecolor='white')
ax1.set_title('$\mu^1$')
# ax1.imshow(mu1[::-1,:], extent=[DomX[0],DomX[-1],DomY[0],DomY[-1]],aspect=DeltaX/DeltaY,vmin=0,vmax=MAX)    
Ratio=mu1/(mu0+mu1)
Ratio=np.nan_to_num(Ratio)
plt.imshow(Ratio[::-1,:], extent=[DomX[0],DomX[-1],DomY[0],DomY[-1]],aspect=0.5*DeltaX/DeltaY)    
plt.plot(Trench["x"],Trench["y"],color="hotpink")
plt.xlim(DomX[0],DomX[-1])
plt.ylim(DomY[0],DomY[-1])
plt.colorbar()
plt.show()


Ratio=mu2/(mu0+mu2)
Ratio=np.nan_to_num(Ratio)
ax2=p2.plot(facecolor="none", edgecolor='white')
ax2.set_title('$\mu^2$')
# ax2.imshow(mu2[::-1,:], extent=[DomX[0],DomX[-1],DomY[0],DomY[-1]],aspect=DeltaX/DeltaY,vmin=0,vmax=MAX)
plt.imshow(Ratio[::-1,:], extent=[DomX[0],DomX[-1],DomY[0],DomY[-1]],aspect=0.5*DeltaX/DeltaY)    
plt.plot(Trench["x"],Trench["y"],color="hotpink")
plt.xlim(DomX[0],DomX[-1])
plt.ylim(DomY[0],DomY[-1])
plt.ylabel("Latitude")
plt.xlabel("Longitude")
plt.colorbar()
plt.show()


Ratio=mu3/(mu0+mu3)
Ratio=np.nan_to_num(Ratio)
ax3=p3.plot(facecolor="none", edgecolor='white')
ax3.set_title('$\mu^3$')
# ax3.imshow(mu3[::-1,:], extent=[DomX[0],DomX[-1],DomY[0],DomY[-1]],aspect=DeltaX/DeltaY,vmin=0,vmax=MAX)
plt.imshow(Ratio[::-1,:], extent=[DomX[0],DomX[-1],DomY[0],DomY[-1]],aspect=0.5*DeltaX/DeltaY)    
plt.plot(Trench["x"],Trench["y"],color="hotpink")
plt.xlim(DomX[0],DomX[-1])
plt.ylim(DomY[0],DomY[-1])
plt.xlabel("Longitude")
plt.colorbar()
plt.show()

Ratio=mu4/(mu0+mu4)
Ratio=np.nan_to_num(Ratio)
ax4=p4.plot(facecolor="none", edgecolor='white')
ax4.set_title('$\mu^4$')
# ax4.imshow(mu4[::-1,:], extent=[DomX[0],DomX[-1],DomY[0],DomY[-1]],aspect=DeltaX/DeltaY,vmin=0,vmax=MAX)
plt.imshow(Ratio[::-1,:], extent=[DomX[0],DomX[-1],DomY[0],DomY[-1]],aspect=0.5*DeltaX/DeltaY)    
plt.plot(Trench["x"],Trench["y"],color="hotpink")
plt.xlim(DomX[0],DomX[-1])
plt.ylim(DomY[0],DomY[-1])
plt.xlabel("Longitude")
plt.colorbar()
plt.show()






#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################


Extension='./figures/'
PII=np.loadtxt(Extension+'PII.csv', delimiter=',')
PIJ=np.load(Extension+'/PIJ.npz')
print(PIJ.files)
PIJ=PIJ['PIJ']
mu0=np.loadtxt(Extension+'mu0.csv', delimiter=',')
mu1=np.loadtxt(Extension+'mu1.csv', delimiter=',')
mu2=np.loadtxt(Extension+'mu2.csv', delimiter=',')
mu3=np.loadtxt(Extension+'mu3.csv', delimiter=',')
mu4=np.loadtxt(Extension+'mu4.csv', delimiter=',')


import shapefile as shp

M=Datos["Magnitude"]





import geopandas as gpd
# Trench=pd.read_csv('./trench_cocos.csv', delimiter=',')
# M=Datos["Magnitude"]


Long=Datos["Longitude"][:]
Lat=Datos["Latitude"][:]


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



file_path = './Cocos_Plate_Z120.dat'
df = pd.read_csv(file_path, delimiter='\t')
df=df.loc[df[df.columns[0]]>DomX[0]] 
df=df.loc[df[df.columns[0]]<DomX[-1]]
df=df.loc[df[df.columns[1]]>DomY[0]]
df=df.loc[df[df.columns[1]]<DomY[-1]] 
np.min(df[df.columns[1]])
Xslab, Yslab = np.meshgrid(np.arange(lonmin,lonmax,0.01), np.arange(latmin,latmax,0.01))
Zslab = griddata((df[df.columns[0]], df[df.columns[1]]), df[df.columns[2]], (Xslab, Yslab), method='cubic')
np.max(df[df.columns[2]])



shape=shp.Reader("./Mexico.shp") #http://geoportal.conabio.gob.mx/metadatos/doc/html/dest_2010gw.html

Trench=pd.read_csv('./trench_cocos.csv', delimiter=',')
Reference=np.array([[-104.5,16],[-104,16],[-103.5,16],[-103,16]])
MRef=np.array([4.5,5.5,6.5,7.5])


fig, ax = plt.subplots()

fig.gca().set_facecolor('lightblue')
ax.set_xlim(DomX[0],DomX[-1])
ax.set_ylim(DomY[0],DomY[-1])
for shape in shape.shapeRecords():
    x = [i[0] for i in shape.shape.points[:]]
    y = [i[1] for i in shape.shape.points[:]]
    # plt.plot(x,y,color="gray")
    plt.fill(x, y, color="yellowgreen", edgecolor="lightgray")
    
# gdf = gpd.read_file('/Users/isaias/Desktop/GNSS/Mexico/contdv1mgw.shp')  #https://data.humdata.org/dataset/cod-xa-jpn
# # Set up the figure and axes
# fig, ax = plt.subplots()
# ax.set_facecolor('lightblue')  # Corrected this line
# ax.set_xlim(DomX[0], DomX[-1])
# ax.set_ylim(DomY[0], DomY[-1])

# # Plot the GeoDataFrame on the existing Axes
# gdf.plot(ax=ax, edgecolor='lightgray', facecolor='yellowgreen')
    
    
    


contours=plt.contour(Xslab, Yslab, -Zslab,levels=4,vmin=0, colors='darkviolet')
plt.clabel(contours, inline=True, fontsize=12)
# plt.colorbar()   

gdf = gpd.read_file("FFM2003.geojson")  # Replace with your file path
# gdf.plot(ax=ax,column='slip',edgecolor='black', cmap='hot',figsize=(10, 10))

centroids = gdf.geometry.centroid
x = centroids.x.values
y = centroids.y.values
z = gdf['slip'].values  # Replace with your actual data column
# Create a grid for interpolation
xi = np.linspace(x.min(), x.max(), 200)
yi = np.linspace(y.min(), y.max(), 200)
xi, yi = np.meshgrid(xi, yi)
zi = griddata((x, y), z, (xi, yi), method='cubic')
contours = ax.contour(xi, yi, zi, levels=np.arange(1,6.5,1), cmap='hot',alpha=0.5)
# ax.clabel(contours, inline=True, fontsize=8)
# plt.show()


# 
gdf2 = gpd.read_file("FFM2012.geojson")  # Replace with your file path
centroids = gdf2.geometry.centroid
x = centroids.x.values
y = centroids.y.values
z = gdf2['slip'].values  # Replace with your actual data column
# Create a grid for interpolation
xi = np.linspace(x.min(), x.max(), 200)
yi = np.linspace(y.min(), y.max(), 200)
xi, yi = np.meshgrid(xi, yi)
zi = griddata((x, y), z, (xi, yi), method='cubic')
contours = ax.contour(xi, yi, zi, levels=np.arange(1,6.5,1), cmap='hot',alpha=0.5)
# 
gdf3 = gpd.read_file("FFM2014.geojson")  # Replace with your file path
centroids = gdf3.geometry.centroid
x = centroids.x.values
y = centroids.y.values
z = gdf3['slip'].values  # Replace with your actual data column
# Create a grid for interpolation
xi = np.linspace(x.min(), x.max(), 200)
yi = np.linspace(y.min(), y.max(), 200)
xi, yi = np.meshgrid(xi, yi)
zi = griddata((x, y), z, (xi, yi), method='cubic')
contours = ax.contour(xi, yi, zi, levels=np.arange(1,6.5,1), cmap='hot',alpha=0.5)

##################################################

p=gpd.GeoSeries(df1.geometry[4])
p2=gpd.GeoSeries(df2.geometry[3])
p3=gpd.GeoSeries(df3.geometry[3])
p4=gpd.GeoSeries(df4.geometry)

p.plot(ax=ax,alpha=0.8,facecolor='none',edgecolor="blue")

p2.plot(ax=ax,alpha=0.8,edgecolor="orange",facecolor='none')
p3.plot(ax=ax,alpha=0.8,edgecolor="fuchsia",facecolor='none')
p4.plot(ax=ax,alpha=0.8,edgecolor="black",facecolor='none')
##################################################

PIJReduc=PIJ[:,:][:,:]
for r in range(len(PIJReduc)):
    for s in range(len(PIJReduc)):
        if PIJReduc[s,r]>=0.5:
            ax.annotate("", xy=(Long[Long.index[s]],Lat[Long.index[s]]), xytext=(Long[Long.index[r]],Lat[Long.index[r]]),arrowprops=dict(arrowstyle="->,head_width=0.1",alpha=PIJReduc[s,r]))

    


for i in range(len(DomX)):
    plt.axvline(DomX[i],color="firebrick",alpha=0.5,linewidth=0.5)
    plt.axhline(DomY[i],color="firebrick",alpha=0.5,linewidth=0.5)
plt.plot(Trench["x"],Trench["y"],color="hotpink")
plt.scatter(Datos["Longitude"],Datos["Latitude"],s=np.exp(5*(M-np.min(M))/(np.max(M)-np.min(M))+1),alpha=0.5,color="dodgerblue")
# plt.scatter(Earthquakes2["Longitude"],Earthquakes2["Latitude"],s=np.exp(5*(Earthquakes2["Mw*"]-np.min(Earthquakes2["Mw*"]))/(np.max(Earthquakes2["Mw*"])-np.min(Earthquakes2["Mw*"]))+1),alpha=0.5,color="dodgerblue")


plt.scatter(Reference[:,0],Reference[:,1],s=np.exp(5*(MRef-np.min(M))/(np.max(M)-np.min(M))+1),alpha=0.5,color="Green")
for i in range(len(Reference)):
    plt.text(Reference[i,0]-0.2,Reference[i,1]+0.2, str(i+4.5))

plt.ylabel("Latitude")
plt.xlabel("Longitude")


plt.scatter(-100.2672,17.0485, color="gold",marker="x")

# plt.scatter(Datos.Longitude[27],Datos.Latitude[27], color="gold",marker="x")
# plt.scatter(Datos.Longitude[30],Datos.Latitude[30], color="gold",marker="x")

plt.show()




###################################################################################################
Datos.loc[ (Datos.Longitude<-102) & (Datos.Longitude>-103) & (Datos.Latitude>17.3)& (Datos.Latitude<18.6) & (Datos.Magnitude>5.)]
Datos.loc[[30]]
Datos.loc[28]
Fechas[30]
Fechas[28]
###################################################################################################

#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################





# Extension='./figures/'
# PII=np.loadtxt(Extension+'PII.csv', delimiter=',')
# PIJ=np.load(Extension+'/PIJ.npz')
# print(PIJ.files)
# PIJ=PIJ['PIJ']
# mu0=np.loadtxt(Extension+'mu0.csv', delimiter=',')
# mu1=np.loadtxt(Extension+'mu1.csv', delimiter=',')
# mu2=np.loadtxt(Extension+'mu2.csv', delimiter=',')
# mu3=np.loadtxt(Extension+'mu3.csv', delimiter=',')
# mu4=np.loadtxt(Extension+'mu4.csv', delimiter=',')
# thetaopt=np.loadtxt(Extension+'thetaopt.csv', delimiter=',')

# A,alpha,c,p,d=thetaopt
# print(thetaopt)


# suma=0
# suma=np.sum(mu0*T0[0]+mu1*T0[1]+mu2*T0[2]+mu3*T0[3]+mu4*T0[4])*(DeltaX*DeltaY)
# for i in Datos.index[1:]:
#     xObs,yObs,MObs=Datos.loc[i][["Latitude","Longitude","Magnitude"]]
#     t=Fechas[i]
#     intg=(1-mp.power(c,(p-1))*(((maxdate-t).total_seconds())/sectoday+c)**(1-p))
#     intf=fintegrated(j,alpha,d,xObs,yObs)
#     suma+=klam(MObs,A,alpha,M0)*intg*intf
#     # print(klam(MObs,A,alpha,M0)*intg*intf)

# print(suma,len(PIJ))
# # len(Datos)

# d

# d=0.02
# d=0.002

# d=0.01
# 2*np.sqrt(d)*110
# # 2*np.sqrt(10**-3)*110


# d1=(10**(-2.44+0.59*4.3)/110/2)**2



















    
    
    
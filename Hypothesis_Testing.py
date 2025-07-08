#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  4 13:03:08 2025

@author: isaias
"""



from scipy.stats import invwishart,multivariate_normal,gamma,dirichlet,poisson

SampleSize=100
#########################################Caso Mexico


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


# lammax=np.max((mu0+mu1,mu0+mu2,mu0+mu3,mu0+mu4))


# len(PIJ)

# Eventos=poisson.rvs(T0[0]*lammax*DeltaX/DeltaY,size=1)[0]
# Eventos=poisson.rvs(T0[0]*lammax*(DeltaX*DeltaY),size=1)[0]
# Candidatos=np.vstack((np.random.uniform(Datos["Longitude"].min(),Datos["Longitude"].max(),Eventos),np.random.uniform(Datos["Latitude"].min(),Datos["Latitude"].max(),Eventos),np.sort(np.random.uniform(0,(maxdate-mindate).total_seconds(),Eventos)))).T

# NT=np.zeros(SampleSize)

# def SimulaSintetico(mu0,ind):
#     Simulados=np.zeros(3)
#     if ind==0:
#         Tmin=0
#         Tmax=(maxdate-mindate).total_seconds()
#     else:
#         Tmin=(GNSSDate[DominioInt][2*(ind-1)]-mindate).total_seconds()
#         Tmax=(GNSSDate[DominioInt][2*(ind-1)+1]-mindate).total_seconds()        
    
#     for i in range(len(mu0)):
#         for j in range(len(mu0)):        
#             simular=poisson.rvs(mu0[::-1,:][len(mu0)-1-i,len(mu0)-1-j]*T0[ind]*(DeltaX*DeltaY))
#             if simular>0:
#                 try:
#                     Nuevos=np.vstack((np.random.uniform(DomX[len(mu0)-j],DomX[len(mu0)-j-1],size=simular),np.random.uniform(DomY[i],DomY[i+1],size=simular) ,np.random.uniform(Tmin,Tmax,size=simular))).T
#                     Simulados=np.vstack((Simulados,Nuevos))
#                 except:
#                     pass
    
#     return Simulados[1:,:]

# def SimulaGR(Simulados):
#     return np.reshape(M0-np.log10(1-np.random.uniform(size=len(Simulados)))/bvalue,(len(Simulados),1))

# # def SimulaPL(Simulados):
# #     # return (c/(np.random.uniform(size=len(Simulados)))**(p-1)-c)*sectoday
# #     return c * ((1 - np.random.uniform(size=len(Simulados))) ** (-1 / (p - 1)) - 1)


# def SimulaPL(Simulados):
#     # return (c/(np.random.uniform(size=len(Simulados)))**(p-1)-c)*sectoday
#     return c * ((1 - np.random.uniform(size=len(Simulados))) ** (-1 / (p - 1)) - 1)*sectoday



# for h in range(SampleSize):
    
#     # np.random.seed(1)
#     Simulados=SimulaSintetico(mu0,0)

#     try :
#         Simulados=np.vstack((Simulados,SimulaSintetico(mu1,1)))
#     except:
#         pass        
#     try :
#         Simulados=np.vstack((Simulados,SimulaSintetico(mu2,2)))
#     except:
#         pass  
#     try :
#         Simulados=np.vstack((Simulados,SimulaSintetico(mu3,3)))
#     except:
#         pass  
#     try :
#         Simulados=np.vstack((Simulados,SimulaSintetico(mu4,4)))
#     except:
#         pass  
        
#     # Simulados=np.vstack((SimulaSintetico(mu0,0),SimulaSintetico(mu1,1),SimulaSintetico(mu2,2),SimulaSintetico(mu3,3),SimulaSintetico(mu4,4)))
#     # plt.scatter(Simulados[:,0],Simulados[:,1])
#     # plt.show()
    
#     # np.random.seed(1)
#     Simulados=np.hstack((Simulados,SimulaGR(Simulados), np.zeros((len(Simulados),1))))
#     # plt.hist(Simulados[:,3])
#     # plt.show()
    
#     i=0
#     while i<len(Simulados):
#         if Simulados[i,4]==0:
#             simular=poisson.rvs(np.float32(klam(Simulados[i,3],A,alpha,M0)))
#             Simulados[i,4]=1
#             if simular>0:
#                 Estados=np.zeros((simular,1))
#                 cov=(d*np.exp(alpha*(Simulados[i,3]-M0)))*np.diag((1,1))
#                 Simulados=np.vstack((Simulados,np.hstack((np.reshape(multivariate_normal.rvs(size=simular, mean=(Simulados[i,0],Simulados[i,1]), cov=cov),(simular,2)),np.reshape(Simulados[i,2]+SimulaPL(Estados),(simular,1)),SimulaGR(Estados),Estados))))
#         i+=1
    
#     # plt.scatter(Simulados[:,0],Simulados[:,1])
#     # plt.show()
#     Simulados=Simulados[Simulados[:, 2].argsort()]
#     recuperafecha=np.vectorize(lambda j :mindate+ timedelta(seconds=np.min((Simulados[j,2],(maxdate-mindate).total_seconds()))))

#     recuperafecha=recuperafecha(np.arange(len(Simulados)))

#     ArregloFechas=np.zeros((len(Simulados),6))
#     for j in range(len(Simulados)):
#         ArregloFechas[j,:]=recuperafecha[j].year,recuperafecha[j].month,recuperafecha[j].day,recuperafecha[j].hour,recuperafecha[j].minute,recuperafecha[j].second
    
#     Simulados=np.hstack((Simulados,ArregloFechas))
    
#     # print(len(Datos),len(Simulados))
        
#     Simulados=pd.DataFrame(Simulados)
#     Simulados.columns=["Longitude","Latitude","A","Magnitude","B","Year","Month","Day","Hour","Minute","Second"]
#     Simulados=Simulados.drop(columns=['A', 'B'])
#     Simulados=Simulados.astype({'Year': 'int', "Month": 'int', "Day": 'int',"Hour": 'int', "Minute": 'int',"Second": 'int'})
#     # Simulados=Simulados[Simulados.apply(inside, axis=1)]
#     Simulados=Simulados.loc[ Simulados["Year"]<anio]
    
#     Simulados=Simulados.loc[  (Simulados["Longitude"]>Datos["Longitude"].min()) &   (Simulados["Longitude"]<Datos["Longitude"].max()) ]
#     Simulados=Simulados.loc[  (Simulados["Latitude"]>Datos["Latitude"].min()) &   (Simulados["Latitude"]<Datos["Latitude"].max()) ]
    
#     Simulados=Simulados.reset_index()
#     Simulados.to_csv('./Hypothesis/Simulados'+str(h)+'.csv', index=False)
#     print(len(Simulados),len(Datos))
#     NT[h]=len(Simulados)


# (np.sum((NT-len(Datos))**2),len(Datos)-np.mean(NT))
# with open('./SquareSize.txt', 'a') as file:
#     file.write(str(npartX)+","+str(np.sum((NT-len(Datos))**2))+","+str(len(Datos)-np.mean(NT))+"\n")



# CV=np.loadtxt('./SquareSize.txt',delimiter=",") 
# plt.plot(CV[:,0],CV[:,1]/100)
# plt.show()
####################

###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################


# # for h in np.arange(SampleSize):
# SampleSize=100


# for h in np.arange(66,SampleSize):        
#     Datos=pd.read_csv('./Hypothesis/Simulados'+str(h)+'.csv',delimiter=",")
#     print(len(Datos))
#     bvalue=1/np.mean(Datos["Magnitude"]-M0)*np.log10(np.exp(1))
#     Datos.loc[Datos["Magnitude"]>7]
#     for i in range(len(Datos)):
#         if i==0:
#             Fechas=datetime(Datos["Year"][i], Datos["Month"][i], Datos["Day"][i], Datos["Hour"][i], Datos["Minute"][i], int(Datos["Second"][i]), int(str(Datos["Second"][i])[-1])*100000)
#         else :
#             Fechas=np.hstack((Fechas,datetime(Datos["Year"][i], Datos["Month"][i], Datos["Day"][i], Datos["Hour"][i], Datos["Minute"][i], int(Datos["Second"][i]), int(str(Datos["Second"][i])[-1])*100000)))
#     Diferencias=np.zeros(len(Fechas))
#     for i in range(len(Diferencias)):
#         Diferencias[i]=(Fechas[i]-mindate).total_seconds()    
#     ###############################################################################Algorithm 1 Fox
#     PII=np.zeros((len(Fechas),len(DominioInt)//2+1))
#     PII[:,0]=1/(np.arange(len(Fechas))+1)
#     for i in range(len(DominioInt)//2):
#         PII[(GNSSDate[DominioInt][2*i]<Fechas) * (GNSSDate[DominioInt][2*i+1]>Fechas),0]=PII[(GNSSDate[DominioInt][2*i]<Fechas) * (GNSSDate[DominioInt][2*i+1]>Fechas),0]/2
#         PII[(GNSSDate[DominioInt][2*i]<Fechas) * (GNSSDate[DominioInt][2*i+1]>Fechas),(i+1)]=PII[(GNSSDate[DominioInt][2*i]<Fechas) * (GNSSDate[DominioInt][2*i+1]>Fechas),0]
    
    
#     ###############################################################################
#     ###############################################################################
#     PIJ=np.zeros((len(Fechas),len(Fechas)))
#     for i in range(len(Fechas)):
#         for j in np.arange(i):
#             PIJ[i,j]=1/(i+1)
#     # Trench=pd.read_csv('./trench_cocos.csv', delimiter=',')
#     # Reference=np.array([[-104.5,16],[-104,16],[-103.5,16],[-103,16]])
#     # MRef=np.array([4.5,5.5,6.5,7.5])
    
#     M=Datos["Magnitude"]
    
    
#     mu0=np.zeros((npartY-1,npartX-1))
#     mu1=np.zeros((npartY-1,npartX-1))
#     mu2=np.zeros((npartY-1,npartX-1))
#     mu3=np.zeros((npartY-1,npartX-1))
#     mu4=np.zeros((npartY-1,npartX-1))
    
        
#     ###############################################################################
#     thetaopt=(0.05, #antes 0.04
#               1.5,
#               0.02,
#               1.5,
#               0.004)
    
#     A,alpha,c,p,d=thetaopt
        
    
#     Const=[[0.,1-10**(-10) ],#10**(0)],#A 
#     [0.5, bvalue*np.log(10)],#alpha
#     [10**(-8), 5],#c
#     [1+10**(-10), 2],#p
#     [0+10**(-10),0.02]]# 10**(-2)]]#d antes era -1 con este cambio no necesito penalizar la esperanza 
    
#     sumaold=0
    
#     for s in range(10):
#         for i in np.arange(mu0.shape[0]):
#             for j in np.arange(mu0.shape[1]):
#                 # print(i,j)
#                 # print(len(Datos[(DomX[j]<Datos["Longitude"])*(DomX[j+1]>=Datos["Longitude"])*(DomY[i]<Datos["Latitude"])*(DomY[i+1]>=Datos["Latitude"])]))
#                 Eleccion=PII[(DomX[j]<Datos["Longitude"])*(DomX[j+1]>=Datos["Longitude"])*(DomY[i]<Datos["Latitude"])*(DomY[i+1]>=Datos["Latitude"])]
#                 #print(Eleccion)
#                 mu0[i,j]=np.sum(Eleccion[:,0])/(DeltaX*DeltaY*T0[0])
#                 mu1[i,j]=np.sum(Eleccion[:,1])/(DeltaX*DeltaY*T0[1])
#                 mu2[i,j]=np.sum(Eleccion[:,2])/(DeltaX*DeltaY*T0[2])
#                 mu3[i,j]=np.sum(Eleccion[:,3])/(DeltaX*DeltaY*T0[3])
#                 mu4[i,j]=np.sum(Eleccion[:,4])/(DeltaX*DeltaY*T0[4])  ####Descomentar si es hasta 2017
#         MAX=np.max(np.hstack((mu0,mu1,mu2,mu3,mu4))) #Falta mu4    
#         MIN=np.min(np.hstack((mu0,mu1,mu2,mu3,mu4))) #Falta mu4    
#         # plt.imshow(mu0[:,::-1], extent=[DomX[0],DomX[1],DomY[0],DomY[1]],aspect=DeltaX/DeltaY)
#         # plt.show()
            
#         plt.hist(np.apply_along_axis(sum, 1,PII))
#         plt.show()
            
#         ###############################################################################
#         Optimizacion=minimize(lambda theta: -likelihhod(theta),thetaopt, bounds=Const, method="Nelder-Mead")#,tol=1e-4 )#"Nelder-Mead"  "Powell" "L-BFGS-B" "Nelder-Mead"
#         print(Optimizacion)
#         A,alpha,c,p,d=Optimizacion["x"]
#         # Optimo=thetaIteracion
#         # A,alpha,c,p,d=Optimo
#         # A,alpha,c,p,d=[0.00531568, 0.6095514,  0.02656009 ,1.99870522 ,0.00608718] 
#         thetaopt=(A,alpha,c,p,d)
        
#         paralell= get_context("fork").Pool(cpu_count())
#         results = paralell.map(lam, range(len(PII)))
#         paralell.close()
        
        
#         for i in range(len(PII)):
#             # print(i/len(PII))
#             xObs,yObs,MObs=Datos.loc[i][["Latitude","Longitude","Magnitude"]]
#             # t=Fechas[i]
#             Bloque=RecuperaBloque(yObs,xObs)  #######Aqui podria haber error
#             muxy,nu=results[i]
#             lamObs=muxy+np.sum(nu)
            
#             PII[i,0]=mu0[Bloque]/(lamObs)#*(PII[i,0])
#             PII[i,1]=mu1[Bloque]/(lamObs)*(PII[i,1]!=0)
#             PII[i,2]=mu2[Bloque]/(lamObs)*(PII[i,2]!=0)
#             PII[i,3]=mu3[Bloque]/(lamObs)*(PII[i,3]!=0)
#             PII[i,4]=mu4[Bloque]/(lamObs)*(PII[i,4]!=0)
#             for j in np.arange(i):
#                 PIJ[i,j]=nu[j]/lamObs
                
            
#         np.savetxt("./Hypothesis/figures/PII"+str(h)+".csv", PII, delimiter=",")
#         # np.savetxt("./figures/PIJ.csv", PIJ, delimiter=",")
#         np.savez_compressed("./Hypothesis/figures/PIJ"+str(h), PIJ=PIJ)
#         np.savetxt("./Hypothesis/figures/mu0"+str(h)+".csv", mu0, delimiter=",")
#         np.savetxt("./Hypothesis/figures/mu1"+str(h)+".csv", mu1, delimiter=",")
#         np.savetxt("./Hypothesis/figures/mu2"+str(h)+".csv", mu2, delimiter=",")
#         np.savetxt("./Hypothesis/figures/mu3"+str(h)+".csv", mu3, delimiter=",")
#         np.savetxt("./Hypothesis/figures/mu4"+str(h)+".csv", mu4, delimiter=",")
#         np.savetxt("./Hypothesis/figures/thetaopt"+str(h)+".csv", thetaopt, delimiter=",")
        
        
#         # PIJ=np.load('./Hypothesis/figures/PIJ'+str(h)+'.npz')
#         # print(PIJ.files)
#         # PIJ=PIJ['PIJ']
            
            
#         suma=np.sum(mu0*T0[0]+mu1*T0[1]+mu2*T0[2]+mu3*T0[3]+mu4*T0[4])*(DeltaX*DeltaY)
#         for i in Datos.index[1:]:
#             xObs,yObs,MObs=Datos.loc[i][["Latitude","Longitude","Magnitude"]]
#             t=Fechas[i]
#             intg=(1-mp.power(c,(p-1))*(((maxdate-t).total_seconds())/sectoday+c)**(1-p))
#             intf=fintegrated(j,alpha,d,xObs,yObs)
#             suma+=klam(MObs,A,alpha,M0)*intg*intf
        

#         print(suma,len(PIJ))
#         if np.abs(suma-sumaold)<1:
#             break
#         else :
#             sumaold=np.copy(suma)
    

############################################################
############################################################
############################################################
############################################################
############################################################
############################################################    



# fig2 = px.imshow(PII)
# fig2 = px.imshow(np.hstack((PII,np.zeros((len(PII),len(PII))))))
# fig2.show()
#fig2.write_html("./figures/Genealogy"+".html")



h=0
Datos=pd.read_csv('./Hypothesis/Simulados'+str(h)+'.csv',delimiter=",")
Datos=Datos[Datos.apply(inside, axis=1)]
Datos=Datos.loc[ Datos["Year"]<anio]
print(np.min(Datos["Magnitude"]))
Datos=Datos.loc[ Datos["Magnitude"]>=M0]
Datos=Datos.reset_index()
print(len(Datos))
# bvalue=1/np.mean(Datos["Magnitude"]-M0)*np.log10(np.exp(1))
# avalue=np.log10(len(Datos))+bvalue*M0
Datos.loc[Datos["Magnitude"]>7]
Extension='./Hypothesis/figures/'
PII=np.loadtxt(Extension+'PII'+str(h)+'.csv', delimiter=',')
PIJ=np.load(Extension+'/PIJ'+str(h)+'.npz')
print(PIJ.files)
PIJ=PIJ['PIJ']
mu0=np.loadtxt(Extension+'mu0'+str(h)+'.csv', delimiter=',')
mu1=np.loadtxt(Extension+'mu1'+str(h)+'.csv', delimiter=',')
mu2=np.loadtxt(Extension+'mu2'+str(h)+'.csv', delimiter=',')
mu3=np.loadtxt(Extension+'mu3'+str(h)+'.csv', delimiter=',')
mu4=np.loadtxt(Extension+'mu4'+str(h)+'.csv', delimiter=',')
thetaopt=np.loadtxt(Extension+'thetaopt'+str(h)+'.csv', delimiter=',')

###############################################################################################

import geopandas as gpd

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



h=0
Extension="/Users/isaias/Desktop/GNSS/Hypothesis/Figures/"
mu0=np.loadtxt(Extension+'mu0'+str(h)+'.csv', delimiter=',')

SampleSize=50
mu0B=np.zeros((len(mu0),len(mu0),SampleSize))

def Promedio(mu):
    Extension='./Hypothesis/Figures/'
    for h in range(SampleSize):
        mu0B[:,:,h]=np.loadtxt(Extension+mu+str(h)+'.csv', delimiter=',')
    return mu0B
mu0=Promedio("mu0")
mu0=np.mean(mu0,axis=2)

mu="mu0"
mu=Promedio(mu)
mu0m=np.median(mu,axis=2)
fig, ax = plt.subplots()
im=ax.imshow(mu0m[::-1,:], extent=[DomX[0],DomX[-1],DomY[0],DomY[-1]],aspect=0.5*DeltaX/DeltaY)
np.max(mu0)
fig.colorbar(im, ax=ax,fraction=0.03, pad=0.01)
plt.show()


import shapefile as shp

def Promedio(mu):
    Extension='./Hypothesis/Figures/'
    for h in range(SampleSize):
        muaux=np.loadtxt(Extension+mu+str(h)+'.csv', delimiter=',')
        mu0B[:,:,h]=muaux/(np.loadtxt(Extension+"mu0"+str(h)+'.csv', delimiter=',')+muaux)
        mu0B[np.isnan(mu0B)]=0
    return mu0B


level= 0.1#0.05
cutValue=0.
for mu in ["mu1","mu2","mu3","mu4"]:
    # mu0=Promedio(mu)
    # mu0m=np.mean(mu0,axis=2)
    if mu=="mu1":
        mu0m=mu1
    if mu=="mu2":
        mu0m=mu2
    if mu=="mu3":
        mu0m=mu3
    if mu=="mu4":
        mu0m=mu4
    mu0m=mu0m/(mu0+mu0m)
    mu0m[np.isnan(mu0m)]=0

        
    # mu0q=np.quantile(mu0,level,axis=2)
    mu0q=np.quantile(Promedio(mu),level,axis=2)
    mu0m[mu0q<=cutValue]=0
    
    
    
    shape=shp.Reader("./Mexico.shp") #http://geoportal.conabio.gob.mx/metadatos/doc/html/dest_2010gw.html    
    Trench=pd.read_csv('./trench_cocos.csv', delimiter=',')
    
    fig, ax = plt.subplots()
    
    # fig.gca().set_facecolor('lightblue')
    ax.set_xlim(DomX[0],DomX[-1])
    ax.set_ylim(DomY[0],DomY[-1])
    for shape in shape.shapeRecords():
        x = [i[0] for i in shape.shape.points[:]]
        y = [i[1] for i in shape.shape.points[:]]
        # plt.plot(x,y,color="gray")
        plt.fill(x, y, color="none", edgecolor="gray")
    

    
 
    
    # fig, ax = plt.subplots()
    # im=ax.imshow(mu0m[::-1,:], extent=[DomX[0],DomX[-1],DomY[0],DomY[-1]],aspect=0.5*DeltaX/DeltaY)
    # np.max(mu0)
    # fig.colorbar(im, ax=ax,fraction=0.03, pad=0.01)
    # plt.show()
    im=ax.imshow(mu0m[::-1,:], extent=[DomX[0],DomX[-1],DomY[0],DomY[-1]],aspect=0.5*DeltaX/DeltaY,vmin=0,vmax=1)
    if mu=="mu1":
        p.plot(ax=ax,facecolor="none", edgecolor='white')
    if mu=="mu2":
        p2.plot(ax=ax,facecolor="none", edgecolor='white')
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
    if mu=="mu3":
        p3.plot(ax=ax,facecolor="none", edgecolor='white')
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
    if mu=="mu4":
        p4.plot(ax=ax,facecolor="none", edgecolor='white')
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

    
    fig.colorbar(im, ax=ax,fraction=0.03, pad=0.01)
    ax.set_ylabel("Latitude")
    ax.set_xlabel("Longitude")
    ax.plot(Trench["x"],Trench["y"],color="hotpink")
    ax.set_xlim(DomX[0],DomX[-1])
    ax.set_ylim(DomY[0],DomY[-1])
    #plt.scatter(Datos["Longitude"],Datos["Latitude"],s=np.exp(5*(M-np.min(M))/(np.max(M)-np.min(M))+1),alpha=0.1)
    plt.show()




























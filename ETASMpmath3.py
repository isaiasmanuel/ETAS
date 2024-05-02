
from mpmath import mp
import numpy as np
import pandas as pd
import scipy as sp
import geopy.distance
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import time
from multiprocessing import Pool, cpu_count, get_context, Process

    
################ Read data

Datos=pd.read_csv('/Users/isaias/Desktop/hypo_MASUB.csv',delimiter=",")
anio=2003
Datos=Datos.loc[ Datos["Year"]<anio]
Datos.loc[0]


sectoyear=86400

################ Dates and distances

for i in range(len(Datos)):
  if i==0:
    Fechas=datetime(Datos["Year"][i], Datos["Month"][i], Datos["Day"][i], Datos["Hour"][i], Datos["Minute"][i], int(Datos["Second"][i]), int(str(Datos["Second"][i])[-1])*100000)
  else :
    Fechas=np.hstack((Fechas,datetime(Datos["Year"][i], Datos["Month"][i], Datos["Day"][i], Datos["Hour"][i], Datos["Minute"][i], int(Datos["Second"][i]), int(str(Datos["Second"][i])[-1])*100000)))


Dist=sp.spatial.distance_matrix(Datos[["Longitude","Latitude"]], Datos[["Longitude","Latitude"]], p=2, threshold=1000000)
Dist=Dist**2
plt.matshow(Dist)
plt.colorbar()
plt.show()


mindate=np.min(Fechas)
mindate=datetime(2000, 1,1)
maxdate=datetime(anio, 1,1)

############## Integration region

np.max(Fechas)


Tmax=(maxdate-mindate).total_seconds()/sectoyear
latmin=np.min(Datos["Latitude"])-1
latmax=np.max(Datos["Latitude"])+1
lonmin=np.min(Datos["Longitude"])-1
lonmax=np.max(Datos["Longitude"])+1


print(maxdate-mindate)
print(np.min(Datos["Latitude"]),np.max(Datos["Latitude"]))
print(np.min(Datos["Longitude"]),np.max(Datos["Longitude"]))
############ 
Coord=np.vstack((Datos["Latitude"],Datos["Longitude"]))
Coord=np.transpose(Coord)


##########################################################
Diferencias=np.zeros(len(Fechas))

for i in range(len(Diferencias)):
    Diferencias[i]=(Fechas[i]-mindate).total_seconds()

DiferenciasMax=np.zeros(len(Fechas))

for i in range(len(Diferencias)):
    DiferenciasMax[i]=(maxdate-Fechas[i]).total_seconds()



############# Define useful functions


def klam(M,A,alpha,M0):
    return A*mp.exp(alpha*(M-M0))


def f(x,y,M, d,alpha,M0):
  return 1/(2*np.pi*d*np.exp(alpha*(M-M0)))*mp.exp(-(Dist[x,y])/(2*d*np.exp(alpha*(M-M0))))

def fint(x,y,M, d,alpha,M0):
  return 1/(2*np.pi*d*mp.exp(alpha*(M-M0)))*mp.exp(-(x**2+y**2)/(2*d*np.exp(alpha*(M-M0))))



def g(t,p,c):
  return (p-1)*mp.power(c,p-1)*mp.power(t+c,-p)*(t>0)


################# Define two lambda, one is to integrate


def fintegrated(i,alpha,d):
    cov=(d*mp.exp(alpha*(Datos["Magnitude"][i]-M0)))*np.diag((1,1))
    cum1=sp.stats.multivariate_normal.cdf((latmax,lonmax), mean=Coord[i,:], cov=cov)
    cum2=sp.stats.multivariate_normal.cdf((latmin,lonmax), mean=Coord[i,:], cov=cov)
    cum3=sp.stats.multivariate_normal.cdf((latmax,lonmin), mean=Coord[i,:], cov=cov)
    cum4=sp.stats.multivariate_normal.cdf((latmin,lonmin), mean=Coord[i,:], cov=cov)
    
    return cum1-cum2-cum3+cum4



def fintegrated2(i,d):
    cov=d**2*np.diag((1,1))
    cum1=sp.stats.multivariate_normal.cdf((latmax,lonmax), mean=Coord[i,:], cov=cov)
    cum2=sp.stats.multivariate_normal.cdf((latmin,lonmax), mean=Coord[i,:], cov=cov)
    cum3=sp.stats.multivariate_normal.cdf((latmax,lonmin), mean=Coord[i,:], cov=cov)
    cum4=sp.stats.multivariate_normal.cdf((latmin,lonmin), mean=Coord[i,:], cov=cov)
    
    return cum1-cum2-cum3+cum4





def lamInt(obs,theta):
    
    x=Datos["Latitude"].loc[obs]
    y=Datos["Longitude"].loc[obs]
    (nu,A,alpha,c,p,d)=theta
    
    
    value=nu*ubar[obs]
    t= mindate+ timedelta(seconds=Diferencias[obs])
    
    
    for Record in Datos.loc[Fechas<t].index:
      Mag=Datos["Magnitude"].loc[Record]
      Lat=Datos["Latitude"].loc[Record]
      Lon=Datos["Longitude"].loc[Record]
      #value+=klam(Mag,A,alpha,M0)*g((t-Fechas[Record]).total_seconds(),p,c)*(1/(2*np.pi*d*np.exp(alpha*(Mag-M0)))*np.exp(-((x-Lat)**2+(y-Lon)**2)/(2*d*np.exp(alpha*(Mag-M0)))))
      value+=klam(Mag,A,alpha,M0)*g((t-Fechas[Record]).total_seconds()/sectoyear,p,c)*fint(x-Lat,y-Lon,Mag, d,alpha,M0)
      #print(i,value,klam(Mag,A,alpha,M0)*g((t-Fechas[Record]).total_seconds(),p,c)*fint(x-Lat,y-Lon,Mag, d,alpha,M0))
    
    return value




#############Calculate likelihood



def likelihhod(thetalike):
    nu,A,alpha,c,p,d=np.copy(thetalike)

    lam=lambda obs: lamInt(obs,thetalike)
    klik=lambda Mag:klam(Mag,A,alpha,M0)
    fvaluesint=np.vectorize(lambda i: fintegrated(i,alpha,d))(np.arange(len(Datos)))
    
        
    
    global lam2
    
    def lam2(obs):
        lam=lambda obs: lamInt(obs,thetalike)
        return mp.log(lam(obs))
    
    paralell= get_context("fork").Pool(cpu_count())
    results = paralell.map(lam2, np.arange(len(Datos)))
    paralell.close()

    suma=np.sum(results)

    integrated=muintegrated*nu
    integral= mp.fadd(integrated,np.sum( fvaluesint[1:]*(1-mp.power(c,(p-1))*(DiferenciasMax[1:]/sectoyear+c)**(1-p))*np.vectorize(klik)(Datos["Magnitude"])[1:]))
    
    lik=mp.fadd(suma, - integral)
    
        
    return lik

############# Calculate rho
def rhoij(i,j):
  Mi=Datos["Magnitude"].loc[i]
  ti,tj=Fechas[[i,j]]
  Latx=Datos["Latitude"].loc[i]
  Lonx=Datos["Longitude"].loc[i]
  Lat=Datos["Latitude"].loc[j]
  Lon=Datos["Longitude"].loc[j]

  num=(klam(Mi,A,alpha,M0)*g((tj-ti).total_seconds()/sectoyear,p,c)*fint(Lat-Latx,Lon-Lonx,Mi, d,alpha,M0))
  return num
##############Calculate d_j

def dj(j,rad):
  Eval=np.sum(Dist[:,j]<rad)
  if (Eval>=nnp):
    return (Eval-nnp)+ (rad/(1+rad))
  else :
    return len(Datos)

###############


##########################################################

minradius=0.02
nnp=10
M0=4


##########################################################

djs=np.zeros(len(Datos))

for j in range(len(Datos)):
  Optim=sp.optimize.minimize_scalar(lambda rad: dj(j,rad))
  if Optim["fun"]==len(Datos):
    Optim=sp.optimize.minimize_scalar(lambda rad: dj(j,rad), bounds=(minradius, np.max(Dist) ), method='bounded')
    if Optim["fun"]==len(Datos):
      print("Warning "+str(j))
  djs[j]=Optim["x"]

djs[djs<minradius]=minradius
##########################################################



# Const=[]

# for i in range(6):
#     #Const.append([0.001,1000000000.])
#     Const.append([0+10**(-100),np.infty])
# Const[-2][0]=1+10**(-10)


Const=[[0, 1000], #nu
  [0.01, 100 ],#A
  [0.05, 1.5],#alpha
  [0+10**(-3), 10**5],#c
  [1+10**(-10), 10],#p
  [0+10**(-3), 10**(-1)]]#d

######Inicial

def u(x,y):
    return 1
######

thetaopt=(2,1,1.5,0.1,2,20)
# likelihhod(thetaopt)

thetaopt=(1,0.3,0.01,0.89,1.16,0.002)
nu,A,alpha,c,p,d=thetaopt
muintegrated=mp.mpmathify(Tmax*sp.integrate.dblquad(u, latmin, latmax, lonmin, lonmax)[0])
ubar=np.vectorize(lambda i: u(Datos["Latitude"][i],Datos["Longitude"][i]))(np.arange(len(Datos)))

# thetaopt=(10**200,A,alpha,c,p,d)
print(nu,A,alpha,c,p,d)
likelihhod(thetaopt)
# thetalike=np.copy(thetaopt)





for i in range(10):
    #thetaopt=(1,A,alpha,c,p,d)
    #Optimizacion=sp.optimize.minimize(lambda theta: -likelihhod(theta),thetaopt,bounds=Const, method= "Nelder-Mead" )
    thetaopt=(1,A,alpha,c,p,d)
    # thetaopt=(1,0.3,0.01,0.89,1.16,0.002)
    Optimizacion=sp.optimize.minimize(lambda theta: -likelihhod(theta),thetaopt, bounds=Const, method= "Nelder-Mead" )
    print(Optimizacion)
    nu,A,alpha,c,p,d=Optimizacion["x"]
    thetaopt=Optimizacion["x"]
    # lam=lambda t,x,y: lamInt2(t,x,y,thetaopt)
    
    
    rhojs=np.zeros(len(Datos))
    for j in range(len(Datos)):
      rhoj=lambda i: rhoij(i,j)
      if j!=0:
          rhojs[j]=np.sum(np.vectorize(rhoj)([np.arange(j)])) / lamInt(j,thetaopt)
          if rhojs[j]>1:
              print(rhojs[j])
              rhojs[j]=np.min((rhojs[j],1))
          
      else :
          rhojs[j]=0
    print(rhojs)
    
    plt.figure()
    plt.hist(rhojs)
    plt.xlabel("Frequency")
    plt.title(r"histogram of $\rho_{j}$"+"\n"+str(thetaopt))
    plt.show()
    

    
    
    feature_x = np.linspace(np.min(Datos["Latitude"]), np.max(Datos["Latitude"]), 70)
    feature_y = np.linspace(np.min(Datos["Longitude"]),np.max(Datos["Longitude"]), 70)
    [X, Y] = np.meshgrid(feature_x, feature_y) 
    
    j=0
    def u(x,y):
      suma=0
      for j in range(len(Datos)):
        lat=Datos["Latitude"].loc[j]
        lon=Datos["Longitude"].loc[j]
        #suma+=np.exp(-((x-lat)**2+(y-lon)**2)/(2*djs[j]**2))/(2*djs[j]*np.pi)
        suma+=(1-rhojs[j])*mp.exp(-((x-lat)**2+(y-lon)**2)/(2*djs[j]**2))/(2*np.pi*djs[j]**2)
        if suma<0:
            print(j)
    
    
      return suma/Tmax
    
    
    
    uxy=np.zeros((len(feature_x), len(feature_y)))
    
    for i in np.arange(len(feature_x)):
      for j in np.arange(len(feature_y)):
        uxy[i,j]=u(feature_x[i],feature_y[j])
    

    # nu2=np.copy(nu)
    # u2=sp.interpolate.RectBivariateSpline(feature_x, feature_y, uxy)
    # def u(x,y):
    #     return np.max((u2(x,y)[0][0],0))#*nu2
    
    
    
  
    plt.figure()
    plt.contourf(X, Y, np.log10(uxy))
    # plt.contourf(X, Y, uxy)
    plt.colorbar()
    plt.scatter(Datos["Latitude"],Datos["Longitude"])
    plt.title(r"$\log\hat{\mu}(x,y)$")
    #ax.set_title('Filled Contour Plot')
    #ax.set_xlabel('feature_x')
    #ax.set_ylabel('feature_y')
    plt.show()
    
    
    

    # muintegrated=mp.quad(lambda x,y: Tmax*u(x,y), [latmin, latmax], [lonmin, lonmax])
    muintegrated=np.sum(np.vectorize(lambda i: fintegrated2(i,djs[i]))(np.arange(len(Datos)))*(1-rhojs))##Tmax
    print("integral", muintegrated)
    ubar=np.vectorize(lambda i: u(Datos["Latitude"][i],Datos["Longitude"][i]))(np.arange(len(Datos)))
    # mp.mpmathify(Tmax*sp.integrate.dblquad(u, latmin, latmax, lonmin, lonmax)[0])



#######################


def lamt(t):
    t= mindate+ timedelta(seconds=t)
    
    fvaluesint=np.vectorize(lambda i: fintegrated(i,alpha,d))(np.arange(len(Datos)))
    
    value=muintegrated/Tmax
    for Record in Datos.loc[Fechas<t].index:
        Mag=Datos["Magnitude"].loc[Record]
        value+=klam(Mag,A,alpha,M0)*g((t-Fechas[Record]).total_seconds()/sectoyear,p,c)*fvaluesint[Record]
    return float(value)

sop=np.arange(100)*Tmax/100*sectoyear
evalua=np.vectorize(lamt)(sop)
plt.plot(sop,evalua/(muintegrated/Tmax))
plt.ylabel("$p(t)$")
plt.xlabel("$t$")
